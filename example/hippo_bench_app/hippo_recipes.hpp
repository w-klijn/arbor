#pragma once

#include <cstddef>
#include <memory>
#include <stdexcept>

#include <recipe.hpp>
#include <util/optional.hpp>

#include "../miniapp/morphology_pool.hpp"

// miniapp-specific recipes

namespace hippo {

struct probe_distribution {
    float proportion = 1.f; // what proportion of cells should get probes?
    bool all_segments = true;    // false => soma only
    bool membrane_voltage = true;
    bool membrane_current = true;
};

struct basic_recipe_param {
    // `num_compartments` is the number of compartments to place in each
    // unbranched section of the morphology, A value of zero indicates that
    // the number of compartments should equal the number of piecewise
    // linear segments in the morphology description of that branch.
    //unsigned compartments_per_segment = 1;

    // Total number of synapses on each cell.
    //unsigned num_synapses = 1;

    //std::string synapse_type = "expsyn";
    float min_connection_delay_ms = 20.0;
    float mean_connection_delay_ms = 20.75;
    float syn_weight_per_cell = 0.3;

    arb::morphology_pool morphologies = arb::default_morphology_pool;

    // If true, iterate through morphologies rather than select randomly.
    bool morphology_round_robin = false;

    // If set we are importing the spikes injected in the network from file
    // instead of a single spike at t==0
    arb::util::optional<std::string> input_spike_path;  // Path to file with spikes

    arb::util::optional<std::string> json_populations;
    arb::util::optional<std::string> json_connectome;
};


std::unique_ptr<arb::recipe> make_hippo_recipe(
        basic_recipe_param param,
        probe_distribution pdist = probe_distribution{});

} // namespace arb
