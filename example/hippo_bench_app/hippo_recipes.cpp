#include <cmath>
#include <random>
#include <vector>
#include <utility>

#include <cell.hpp>
#include <dss_cell_description.hpp>
#include <event_generator.hpp>
#include <ipss_cell.hpp>
#include <rss_cell.hpp>
#include <morphology.hpp>
#include <common_types.hpp>
#include <util/debug.hpp>


#include "io.hpp"
#include "hippo_recipes.hpp"

#include "../miniapp/morphology_pool.hpp"

#include "../con_gen/connection_generator.hpp"
#include "../con_gen/con_gen_utils.hpp"

using namespace arb;

namespace hippo {

// TODO: split cell description into separate morphology, stimulus, mechanisms etc.
// description for greater data reuse.

template <typename RNG>
cell make_basic_cell(
    const morphology& morph,
    unsigned num_synapses,
    RNG& rng,
    nlohmann::json const& opts)
{
    arb::cell cell = make_cell(morph, true);

    // TODOW: It all fails of the json is not completely specified.
    // It would be nice to have some 'standard' way to handle this

    unsigned compartments_per_segment = opts["compartments_per_segment"];
    std::string syn_type(opts["synapse_type"].get<std::string>());

    for (auto& segment: cell.segments()) {
        if (compartments_per_segment !=0) {
            if (cable_segment* cable = segment->as_cable()) {
                cable->set_compartments(compartments_per_segment);
            }
        }

        if (segment->is_dendrite()) {
            segment->add_mechanism(opts["dendrite_mechanism"].get<std::string>());
            segment->rL = opts["dendrite_rL"];
        }
    }

    cell.soma()->add_mechanism(opts["soma_mechanism"].get<std::string>());
    cell.add_detector({0,0}, 20);

    auto distribution = std::uniform_real_distribution<float>(0.f, 1.0f);

    // TODOW: We might have to add a better synapse location here at some time
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

    arb::mechanism_spec syn_default(syn_type);
    for (unsigned i=0; i<num_synapses; ++i) {
        unsigned id = terminals[i%terminals.size()];


        cell.add_synapse({id, distribution(rng)}, syn_default);
    }

    return cell;
}

class hippo_recipe: public recipe {
public:
    hippo_recipe(hippo::hippo_recipe_param param, hippo::probe_distribution pdist):
        // TODO fix ncless
        param_(std::move(param)), pdist_(std::move(pdist))
    {
        // Cells are not allowed to connect to themselves; hence there must be least two cells
        // to build a connected network.

        EXPECTS(param_.morphologies.size()>0);

        if (param_.json_connectome && param_.json_populations) {
            con_gen_ = arb_con_gen::connection_generator(
                con_gen_util::parse_populations_from_path(param_.json_populations.value()),
                con_gen_util::parse_projections_from_path(param_.json_connectome.value()));
        }
        else {
            con_gen_ = arb_con_gen::connection_generator(
                con_gen_util::default_populations(),
                con_gen_util::default_connectome());
        }

    }

    cell_size_type num_cells() const override {
        return con_gen_.num_cells();
    }

    util::unique_any get_cell_description(cell_gid_type i) const override {
        auto kind = con_gen_.get_cell_kind(i);
        nlohmann::json const& opts = con_gen_.get_cell_opts(i);

        if (kind == arb::cell_kind::inhomogeneous_poisson_spike_source) {

            ipss_cell_description cell(opts);

            return util::unique_any(std::move(cell));
        }
        else if (kind == arb::cell_kind::data_spike_source) {

            std::vector<time_type> spike_times;
            for (auto & entry : opts["spike_times"]) {
                spike_times.push_back({ entry });
            }
            dss_cell_description cell(spike_times);

            return util::unique_any(std::move(cell));
        }
        else if (kind == arb::cell_kind::regular_spike_source) {

            rss_cell cell(opts["start_time"], opts["period"], opts["stop_time"]);

            return util::unique_any(std::move(cell));
        }


        //else (kind == arb::cell_kind::cable1d_neuron
        auto gen = std::mt19937(i); // TODO: replace this with hashing generator...
        const auto& morph = get_morphology(i);
        unsigned cell_segments = morph.components();

        auto cell = make_basic_cell(morph, con_gen_.num_synapses_on(i), gen, opts);

EXPECTS(cell.num_segments() == cell_segments);

return util::unique_any(std::move(cell));

    }

    std::vector<cell_connection> connections_on(cell_gid_type i) const override {
        std::vector<cell_connection> conns;

        auto connections = con_gen_.synapses_on(i);
        unsigned tag_idx = 0;
        for (auto& syn_par : connections) {
            conns.push_back(cell_connection(
            { syn_par.gid, 0 }, { i, tag_idx }, syn_par.weight, syn_par.delay));
            tag_idx++;
        }
        return conns;
    }

    probe_info get_probe(cell_member_type probe_id) const override {
        if (probe_id.index >= num_probes(probe_id.gid)) {
            throw invalid_recipe_error("invalid probe id");
        }

        // if we have both voltage and current probes, then order them
        // voltage compartment 0, current compartment 0, voltage compartment 1, ...

        cell_probe_address::probe_kind kind;

        int stride = pdist_.membrane_voltage + pdist_.membrane_current;

        if (stride == 1) {
            // Just one kind of probe.
            kind = pdist_.membrane_voltage ?
                cell_probe_address::membrane_voltage : cell_probe_address::membrane_current;
        }
        else {
            EXPECTS(stride == 2);
            // Both kinds available.
            kind = (probe_id.index%stride == 0) ?
                cell_probe_address::membrane_voltage : cell_probe_address::membrane_current;
        }

        cell_lid_type compartment = probe_id.index / stride;
        segment_location loc{ compartment, compartment ? 0.5 : 0.0 };

        // Use probe kind as the token to be passed to a sampler.
        return { probe_id, (int)kind, cell_probe_address{loc, kind} };
    }

    cell_kind get_cell_kind(cell_gid_type i) const override {

        return con_gen_.get_cell_kind(i);
    }

    cell_size_type num_sources(cell_gid_type i) const override {
        return 1;
    }

    cell_size_type num_targets(cell_gid_type i) const override {
        nlohmann::json const& opts = con_gen_.get_cell_opts(i);
        return opts["synapses_per_cell"];
    }

    cell_size_type num_probes(cell_gid_type i) const override {
        bool has_probe = (std::floor(i*pdist_.proportion) != std::floor((i - 1.0)*pdist_.proportion));

        if (!has_probe) {
            return 0;
        }
        else {
            cell_size_type np = pdist_.all_segments ? get_morphology(i).components() : 1;
            return np*(pdist_.membrane_voltage + pdist_.membrane_current);
        }
    }

    // Return two generators attached to the one cell.
    std::vector<arb::event_generator_ptr> event_generators(cell_gid_type gid) const override {
        auto kind = con_gen_.get_cell_kind(gid);

        // For now only cable1d_neuron support events
        if (kind != arb::cell_kind::cable1d_neuron) {
            return {};
        }

        nlohmann::json const& opts = con_gen_.get_cell_opts(gid);

        // If there are no "poisson_generators" event generators
        if (opts.find("poisson_generators") == opts.end()) {
            return {};
        }

        // We need a random number generator
        using RNG = std::mt19937_64;
        using pgen = arb::poisson_generator<RNG>;

        auto hz_to_freq = [](double hz) { return hz*1e-3; };


        // a counter needed to assure unique rng seeds per event generator
        unsigned idx = 1;

        std::vector<arb::event_generator_ptr> gens;
        for (nlohmann::json::const_iterator it = opts["poisson_generators"].begin();
            it != opts["poisson_generators"].end(); ++it) {

            // Get values
            double rate = it.value()["rate"];
            double weight = it.value()["weight"];
            double start_time = it.value()["start"];

            // Sanity skip if generator would not do anything
            if (rate == 0.0 || weight == 0.0) {
                continue;
            }

            gens.push_back(
                arb::make_event_generator<pgen>(
                    cell_member_type{ 0,0 },        // Target synapse (gid, local_id).
                    weight,
                    RNG(gid + idx * 29562872),   // TODO: Better random number generation
                    start_time,
                    hz_to_freq(rate)));
            idx++;

        }
        return gens;
    }

protected:

    hippo::hippo_recipe_param param_;

    hippo::probe_distribution pdist_;

    using exp_param = std::exponential_distribution<float>::param_type;

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


std::unique_ptr<arb::recipe> make_hippo_recipe(
    hippo::hippo_recipe_param param,
    hippo::probe_distribution pdist) {
    return std::unique_ptr<recipe>(new hippo_recipe(param, pdist));
}



} // namespace hippo

