#include <cmath>
#include <random>
#include <vector>
#include <utility>

#include <cell.hpp>
#include <dss_cell_description.hpp>
#include <event_generator.hpp>
#include <rss_cell.hpp>
#include <morphology.hpp>
#include <util/debug.hpp>

#include "io.hpp"
#include "hippo_bench_recipes.hpp"
#include "morphology_pool.hpp"
#include "../con_gen/connection_generator.hpp"
#include "../con_gen/con_gen_utils.hpp"



namespace arb {

// TODO: split cell description into separate morphology, stimulus, mechanisms etc.
// description for greater data reuse.

template <typename RNG>
cell make_basic_cell(
    const morphology& morph,
    arb_con_gen::cell_pars const & cell_pars,
    RNG& rng)
{
    arb::cell cell = make_cell(morph, true);

    for (auto& segment: cell.segments()) {
        if (cell_pars.compartments_per_segment!=0) {
            if (cable_segment* cable = segment->as_cable()) {
                cable->set_compartments(cell_pars.compartments_per_segment);
            }
        }

        if (segment->is_dendrite()) {
            segment->add_mechanism(cell_pars.dendrite_mechanism);
            segment->rL = cell_pars.dendrite_rL;
        }
    }

    cell.soma()->add_mechanism(cell_pars.soma_mechanism);
    cell.add_detector({0,0}, cell_pars.spike_threshold);

    auto distribution = std::uniform_real_distribution<float>(0.f, 1.0f);

    // Distribute the synapses at random locations the terminal dendrites in a
    // round robin manner.
    morph.assert_valid();
    std::vector<unsigned> terminals;
    for (const auto& section: morph.sections) {
        // Note that morphology section ids should match up exactly with cell
        // segment ids!
        if (section.terminal) {
            terminals.push_back(section.id);
        }
    }

    EXPECTS(!terminals.empty());

    arb::mechanism_spec syn_default(cell_pars.synapse_type);
    for (unsigned i=0; i<cell_pars.synapses_per_cell; ++i) {
        unsigned id = terminals[i%terminals.size()];
        cell.add_synapse({id, distribution(rng)}, syn_default);
    }

    return cell;
}

class hippo_bench_recipe: public recipe {
public:
    hippo_bench_recipe(cell_gid_type ncell, basic_recipe_param param,
        probe_distribution pdist = probe_distribution{}):
        ncell_(ncell), param_(std::move(param)), pdist_(std::move(pdist))
    {
        // Cells are not allowed to connect to themselves; hence there must be least two cells
        // to build a connected network.
        if (ncell<2) {
            throw std::runtime_error("A randomly connected network must have at least 2 cells.");
        }
        con_gen_ = arb_con_gen::connection_generator(con_gen_util::parse_populations_from_path(
            "../populations.json"), con_gen_util::parse_projections_from_path("../projections.json")
        );


        EXPECTS(param_.morphologies.size()>0);
        delay_distribution_param_ = exp_param{param_.mean_connection_delay_ms
                            - param_.min_connection_delay_ms};
    }

    cell_size_type num_cells() const override {
        return ncell_ + 1;  // We automatically add a fake cell to each recipe!
    }

    util::unique_any get_cell_description(cell_gid_type i) const override {
        // The last 'cell' is a spike source cell. Either a regular spiking
        // or a spikes from file.
        auto gen = std::mt19937(i); // TODO: replace this with hashing generator...

        const auto& morph = get_morphology(i);
        unsigned cell_segments = morph.components();

        auto cell = make_basic_cell(morph, con_gen_.get_cell_opts(i), gen);

        EXPECTS(cell.num_segments()==cell_segments);
        EXPECTS(cell.synapses().size()==num_targets(i));
        EXPECTS(cell.detectors().size()==num_sources(i));

        return util::unique_any(std::move(cell));
    }

    probe_info get_probe(cell_member_type probe_id) const override {
        if (probe_id.index>=num_probes(probe_id.gid)) {
            throw invalid_recipe_error("invalid probe id");
        }

        // if we have both voltage and current probes, then order them
        // voltage compartment 0, current compartment 0, voltage compartment 1, ...

        cell_probe_address::probe_kind kind;

        int stride = pdist_.membrane_voltage+pdist_.membrane_current;

        if (stride==1) {
            // Just one kind of probe.
            kind = pdist_.membrane_voltage?
                cell_probe_address::membrane_voltage: cell_probe_address::membrane_current;
        }
        else {
            EXPECTS(stride==2);
            // Both kinds available.
            kind = (probe_id.index%stride==0)?
                cell_probe_address::membrane_voltage: cell_probe_address::membrane_current;
        }

        cell_lid_type compartment = probe_id.index/stride;
        segment_location loc{compartment, compartment? 0.5: 0.0};

        // Use probe kind as the token to be passed to a sampler.
        return {probe_id, (int)kind, cell_probe_address{loc, kind}};
    }

    cell_kind get_cell_kind(cell_gid_type i) const override {
        return con_gen_.get_cell_kind(i);

    }

    cell_size_type num_sources(cell_gid_type i) const override {
        return 1;
    }

    cell_size_type num_targets(cell_gid_type i) const override {
        return param_.num_synapses;
    }

    cell_size_type num_probes(cell_gid_type i) const override {
        bool has_probe = (std::floor(i*pdist_.proportion)!=std::floor((i-1.0)*pdist_.proportion));

        if (!has_probe) {
            return 0;
        }
        else {
            cell_size_type np = pdist_.all_segments? get_morphology(i).components(): 1;
            return np*(pdist_.membrane_voltage+pdist_.membrane_current);
        }
    }

    // Return two generators attached to the one cell.
    std::vector<arb::event_generator_ptr> event_generators(cell_gid_type gid) const override {
        if (gid == ncell_) {
            return {};
        }
        using RNG = std::mt19937_64;
        using pgen = arb::poisson_generator<RNG>;

        auto hz_to_freq = [](double hz) { return hz*1e-3; };
        time_type t0 = 0;

        std::vector<arb_con_gen::poisson_event_pars> generators;

        generators.push_back({ 500, 0.001, 0 });
        generators.push_back({ 20, -0.005, 0 });

        //generators.push_back({ 500, 0.001, 0 });
        //generators.push_back({ 20, -0.005, 0 });

        //std::vector<arb::event_generator_ptr> gens;
        //if (gid == 0) {
        //    for (auto const& pars : generators) {
        //        std::cout << pars.rate << "," << pars.weight << ", " << pars.start << "\n";


        //    }
        //}
        //    gens.push_back(
        //        arb::make_event_generator<pgen>(
        //            cell_member_type{ gid,0 }, // Target synapse (gid, local_id).
        //            hz_to_freq(pars.rate),                   // Weight of events to deliver
        //            RNG(29562872),         // Random number generator to use
        //            pars.start,                    // Events start being delivered from this time
        //            pars.weight));            // Expected frequency (events per ms)


        //}

        // Make two event generators.
        std::vector<arb::event_generator_ptr> gens;

        // Add excitatory generator
        gens.push_back(
            arb::make_event_generator<pgen>(
                cell_member_type{ 0,0 }, // Target synapse (gid, local_id).
                generators[0].weight,                   // Weight of events to deliver
                RNG(29562872),         // Random number generator to use
                t0,                    // Events start being delivered from this time
                hz_to_freq(generators[0].rate)));            // Expected frequency (events per ms)

                                       // Add inhibitory generator
        gens.push_back(
            arb::make_event_generator<pgen>(
                cell_member_type{ 0,0 }, generators[1].weight, RNG(86543891), t0, hz_to_freq(generators[1].rate)));

        return gens;
    }

    std::vector<cell_connection> connections_on(cell_gid_type i) const override {
        std::vector<cell_connection> conns;

        //// The rss_cell does not have inputs
        //if (i == ncell_) {
        //    return conns;
        //}
        //auto conn_param_gen = std::mt19937(i); // TODO: replace this with hashing generator...
        //auto source_gen = std::mt19937(i*123+457); // ditto

        //std::uniform_int_distribution<cell_gid_type> source_distribution(0, ncell_-2);

        //for (unsigned t=0; t<param_.num_synapses; ++t) {
        //    auto source = source_distribution(source_gen);
        //    if (source>=i) ++source;

        //    cell_connection cc = draw_connection_params(conn_param_gen);
        //    cc.source = {source, 0};
        //    cc.dest = {i, t};
        //    conns.push_back(cc);

        //    // The rss_cell spikes at t=0, with these connections it looks like
        //    // (source % 20) == 0 spikes at that moment.
        //    if (source % 20 == 0) {
        //        cc.source = {ncell_, 0};
        //        conns.push_back(cc);
        //    }
        //}

        return conns;
    }

protected:
    template <typename RNG>
    cell_connection draw_connection_params(RNG& rng) const {
        std::exponential_distribution<float> delay_dist(delay_distribution_param_);
        float delay = param_.min_connection_delay_ms + delay_dist(rng);
        float weight = param_.syn_weight_per_cell/param_.num_synapses;
        return cell_connection{{0, 0}, {0, 0}, weight, delay};
    }

    cell_gid_type ncell_;
    basic_recipe_param param_;
    probe_distribution pdist_;

    using exp_param = std::exponential_distribution<float>::param_type;
    exp_param delay_distribution_param_;

    arb_con_gen::connection_generator con_gen_;

    const morphology& get_morphology(cell_gid_type gid) const {
        // Allocate to gids sequentially?
        if (param_.morphology_round_robin) {
            return param_.morphologies[gid%param_.morphologies.size()];
        }

        // Morphologies are otherwise selected deterministically pseudo-randomly from pool.
        std::uniform_int_distribution<unsigned> morph_select_dist_(0, param_.morphologies.size()-1);

        // TODO: definitely replace this with a random hash!
        auto gen = std::mt19937(gid+0xbad0cafe);
        return param_.morphologies[morph_select_dist_(gen)];
    }
};



std::unique_ptr<recipe> make_hippo_bench_recipe(
        cell_gid_type ncell,
        basic_recipe_param param,
        probe_distribution pdist)
{
    return std::unique_ptr<recipe>(new hippo_bench_recipe(ncell, param, pdist));
}

} // namespace arb