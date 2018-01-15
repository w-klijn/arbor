#include <cmath>
#include <random>
#include <vector>
#include <utility>

#include <cell.hpp>
#include <dss_cell_description.hpp>
#include <event_generator.hpp>
#include <rss_cell.hpp>
#include <morphology.hpp>
#include <common_types.hpp>
#include <util/debug.hpp>


#include "io.hpp"
#include "hippo_recipes.hpp"
#include "../miniapp/morphology_pool.hpp"

#include "../con_gen/connection_generator.hpp"
#include "../con_gen/con_gen_utils.hpp"


namespace arb {

// TODO: split cell description into separate morphology, stimulus, mechanisms etc.
// description for greater data reuse.

template <typename RNG>
cell make_basic_cell(
    const morphology& morph,
    unsigned compartments_per_segment,
    unsigned num_synapses,
    const std::string& syn_type,
    RNG& rng)
{
    arb::cell cell = make_cell(morph, true);

    for (auto& segment: cell.segments()) {
        if (compartments_per_segment!=0) {
            if (cable_segment* cable = segment->as_cable()) {
                cable->set_compartments(compartments_per_segment);
            }
        }

        if (segment->is_dendrite()) {
            segment->add_mechanism("pas");
            segment->rL = 100;
        }
    }

    cell.soma()->add_mechanism("hh");
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
    hippo_recipe(basic_recipe_param param, probe_distribution pdist):
        // TODO fix ncless
        param_(std::move(param)), pdist_(std::move(pdist)),
        con_gen_(con_gen_util::default_populations(), con_gen_util::default_connectome())
    {
        // Cells are not allowed to connect to themselves; hence there must be least two cells
        // to build a connected network.

        EXPECTS(param_.morphologies.size()>0);
        delay_distribution_param_ = exp_param{param_.mean_connection_delay_ms
                            - param_.min_connection_delay_ms};
    }

    cell_size_type num_cells() const override {
        return con_gen_.num_cells();
    }

    util::unique_any get_cell_description(cell_gid_type i) const override {

        auto kind = con_gen_.get_cell_kind(i);

        if (kind == arb::cell_kind::cable1d_neuron) {
            auto gen = std::mt19937(i); // TODO: replace this with hashing generator...
            const auto& morph = get_morphology(i);
            unsigned cell_segments = morph.components();
            auto cell = make_basic_cell(morph, param_.num_compartments, con_gen_.num_synapses_on(i),
                param_.synapse_type, gen);

            EXPECTS(cell.num_segments() == cell_segments);

            return util::unique_any(std::move(cell));
        }
        if (kind == arb::cell_kind::inhomogeneous_poisson_spike_source) {

        }

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

        return cell_kind::cable1d_neuron;
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

    std::vector<event_generator_ptr> event_generators(cell_gid_type) const override {
        return {};
    }

protected:

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


std::unique_ptr<recipe> make_hippo_recipe(
        basic_recipe_param param,
        probe_distribution pdist)
{
    return std::unique_ptr<recipe>(new hippo_recipe(param, pdist));
}



} // namespace arb