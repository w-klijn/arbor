#pragma once

#include <cstdint>
#include <functional>
#include <iterator>
#include <unordered_map>
#include <vector>

#include <algorithms.hpp>
#include <backends/event.hpp>
#include <cell_group.hpp>
#include <cell.hpp>
#include <common_types.hpp>
#include <event_binner.hpp>
#include <event_queue.hpp>
#include <recipe.hpp>
#include <sampler_map.hpp>
#include <sampling.hpp>
#include <spike.hpp>
#include <util/debug.hpp>
#include <util/double_buffer.hpp>
#include <util/filter.hpp>
#include <util/partition.hpp>
#include <util/range.hpp>
#include <util/unique_any.hpp>

#include <profiling/profiler.hpp>

namespace arb {

template <typename LoweredCell>
class mc_cell_group: public cell_group {
public:
    using lowered_cell_type = LoweredCell;
    using value_type = typename lowered_cell_type::value_type;
    using size_type  = typename lowered_cell_type::value_type;

    mc_cell_group() = default;

    mc_cell_group(std::vector<cell_gid_type> gids, const recipe& rec):
        gids_(std::move(gids))
    {
        // Default to no binning of events
        set_binning_policy(binning_kind::none, 0);

        // Build lookup table for gid to local index.
        for (auto i: util::make_span(0, gids_.size())) {
            gid_index_map_[gids_[i]] = i;
        }

        // Create lookup structure for target ids.
        build_target_handle_partition(rec);
        std::size_t n_targets = target_handle_divisions_.back();

        // Pre-allocate space to store handles, probe map.
        auto n_probes = util::sum_by(gids_, [&rec](cell_gid_type i) { return rec.num_probes(i); });
        probe_map_.reserve(n_probes);
        target_handles_.reserve(n_targets);

        // Construct cell implementation, retrieving handles and maps.
        lowered_.initialize(gids_, rec, target_handles_, probe_map_);

        // Create a list of the global identifiers for the spike sources
        for (auto source_gid: gids_) {
            for (cell_lid_type lid = 0; lid<rec.num_sources(source_gid); ++lid) {
                spike_sources_.push_back({source_gid, lid});
            }
        }
        spike_sources_.shrink_to_fit();
    }

    cell_kind get_cell_kind() const override {
        return cell_kind::cable1d_neuron;
    }

    void reset() override {
        spikes_.clear();
        reset_samplers();
        for (auto& b: binners_) {
            b.reset();
        }
        lowered_.reset();
    }

    void set_binning_policy(binning_kind policy, time_type bin_interval) override {
        binners_.clear();
        binners_.resize(gids_.size(), event_binner(policy, bin_interval));
    }

    void advance(epoch ep, time_type dt, const event_lane_subrange& event_lanes) override {
        PE("advance");
        EXPECTS(lowered_.state_synchronized());
        time_type tstart = lowered_.min_time();

        PE("event-setup");
        staged_events_.clear();
        // skip event binning if empty lanes are passed
        if (event_lanes.size()) {
            for (auto lid: util::make_span(0, gids_.size())) {
                auto& lane = event_lanes[lid];
                for (auto e: lane) {
                    if (e.time>=ep.tfinal) break;
                    e.time = binners_[lid].bin(e.time, tstart);
                    auto h = target_handles_[target_handle_divisions_[lid]+e.target.index];
                    auto ev = deliverable_event(e.time, h, e.weight);
                    staged_events_.push_back(ev);
                }
            }
        }
        PL();

        // Create sample events and delivery information.
        //
        // For each (schedule, sampler, probe set) in the sampler association
        // map that will be triggered in this integration interval, create
        // sample events for the lowered cell, one for each scheduled sample
        // time and probe in the probe set.
        //
        // Each event is associated with an offset into the sample data and
        // time buffers; these are assigned contiguously such that one call to
        // a sampler callback can be represented by a `sampler_call_info`
        // value as defined below, grouping together all the samples of the
        // same probe for this callback in this association.

        struct sampler_call_info {
            sampler_function sampler;
            cell_member_type probe_id;
            probe_tag tag;

            // Offsets are into lowered cell sample time and event arrays.
            sample_size_type begin_offset;
            sample_size_type end_offset;
        };

        PE("sample-event-setup");
        std::vector<sampler_call_info> call_info;

        std::vector<sample_event> sample_events;
        sample_size_type n_samples = 0;
        sample_size_type max_samples_per_call = 0;

        for (auto& sa: sampler_map_) {
            auto sample_times = sa.sched.events(tstart, ep.tfinal);
            if (sample_times.empty()) {
                continue;
            }

            sample_size_type n_times = sample_times.size();
            max_samples_per_call = std::max(max_samples_per_call, n_times);

            for (cell_member_type pid: sa.probe_ids) {
                auto cell_index = gid_to_index(pid.gid);
                auto p = probe_map_[pid];

                call_info.push_back({sa.sampler, pid, p.tag, n_samples, n_samples+n_times});

                for (auto t: sample_times) {
                    sample_event ev{t, cell_index, {p.handle, n_samples++}};
                    sample_events.push_back(ev);
                }
            }
        }

        // Sample events must be ordered by time for the lowered cell.
        util::sort_by(sample_events, [](const sample_event& ev) { return event_time(ev); });
        PL();

        // Run integration.
        lowered_.setup_integration(ep.tfinal, dt, staged_events_, std::move(sample_events));
        PE("integrator-steps");

        while (!lowered_.integration_complete()) {
            lowered_.step_integration();
            if (util::is_debug_mode() && !lowered_.is_physical_solution()) {
                std::cerr << "warning: solution out of bounds  at (max) t "
                          << lowered_.max_time() << " ms\n";
            }
        }
        PL();


        // For each sampler callback registered in `call_info`, construct the
        // vector of sample entries from the lowered cell sample times and values
        // and then call the callback.

        PE("sample-deliver");
        std::vector<sample_record> sample_records;
        sample_records.reserve(max_samples_per_call);

        auto sample_time = lowered_.sample_time();
        auto sample_value = lowered_.sample_value();

        for (auto& sc: call_info) {
            sample_records.clear();
            for (auto i = sc.begin_offset; i!=sc.end_offset; ++i) {
                sample_records.push_back(sample_record{time_type(sample_time[i]), &sample_value[i]});
            }

            sc.sampler(sc.probe_id, sc.tag, sc.end_offset-sc.begin_offset, sample_records.data());
        }
        PL();

        // Copy out spike voltage threshold crossings from the back end, then
        // generate spikes with global spike source ids. The threshold crossings
        // record the local spike source index, which must be converted to a
        // global index for spike communication.

        PE("spike-retrieve");
        for (auto c: lowered_.get_spikes()) {
            spikes_.push_back({spike_sources_[c.index], time_type(c.time)});
        }

        // Now that the spikes have been generated, clear the old crossings
        // to get ready to record spikes from the next integration period.

        lowered_.clear_spikes();
        PL();

        PL();
    }

    // TODOW: Quick hack to be able to set the starting membrane voltage
    void set_resting_potential(double v) {
        lowered_.resting_potential(v);
        lowered_.reset();
    }

    const std::vector<spike>& spikes() const override {
        return spikes_;
    }

    void clear_spikes() override {
        spikes_.clear();
    }

    const std::vector<cell_member_type>& spike_sources() const {
        return spike_sources_;
    }

    void add_sampler(sampler_association_handle h, cell_member_predicate probe_ids,
                     schedule sched, sampler_function fn, sampling_policy policy) override
    {
        std::vector<cell_member_type> probeset =
            util::assign_from(util::filter(util::keys(probe_map_), probe_ids));

        if (!probeset.empty()) {
            sampler_map_.add(h, sampler_association{std::move(sched), std::move(fn), std::move(probeset)});
        }
    }

    void remove_sampler(sampler_association_handle h) override {
        sampler_map_.remove(h);
    }

    void remove_all_samplers() override {
        sampler_map_.clear();
    }

private:
    // List of the gids of the cells in the group.
    std::vector<cell_gid_type> gids_;

    // Hash table for converting gid to local index
    std::unordered_map<cell_gid_type, cell_gid_type> gid_index_map_;

    // The lowered cell state (e.g. FVM) of the cell.
    lowered_cell_type lowered_;

    // Spike detectors attached to the cell.
    std::vector<cell_member_type> spike_sources_;

    // Spikes that are generated.
    std::vector<spike> spikes_;

    // Event time binning manager.
    std::vector<event_binner> binners_;

    // List of events to deliver
    std::vector<deliverable_event> staged_events_;

    // Pending samples to be taken.
    event_queue<sample_event> sample_events_;

    // Handles for accessing lowered cell.
    using target_handle = typename lowered_cell_type::target_handle;
    std::vector<target_handle> target_handles_;

    // Maps probe ids to probe handles (from lowered cell) and tags (from probe descriptions).
    using probe_handle = typename lowered_cell_type::probe_handle;
    probe_association_map<probe_handle> probe_map_;

    // Collection of samplers to be run against probes in this group.
    sampler_association_map sampler_map_;

    // Lookup table for target ids -> local target handle indices.
    std::vector<std::size_t> target_handle_divisions_;

    // Build handle index lookup tables.
    void build_target_handle_partition(const recipe& rec) {
        util::make_partition(target_handle_divisions_,
            util::transform_view(gids_, [&rec](cell_gid_type i) { return rec.num_targets(i); }));
    }

    // Get target handle from target id.
    target_handle get_target_handle(cell_member_type id) const {
        return target_handles_[target_handle_divisions_[gid_to_index(id.gid)]+id.index];
    }

    void reset_samplers() {
        // clear all pending sample events and reset to start at time 0
        sample_events_.clear();
        for (auto &assoc: sampler_map_) {
            assoc.sched.reset();
        }
    }

    cell_gid_type gid_to_index(cell_gid_type gid) const {
        auto it = gid_index_map_.find(gid);
        EXPECTS(it!=gid_index_map_.end());
        return it->second;
    }
};

} // namespace arb
