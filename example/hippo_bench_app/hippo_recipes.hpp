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

struct hippo_recipe_param {

    //std::string synapse_type = "expsyn";
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
        hippo_recipe_param param,
        probe_distribution pdist = probe_distribution{});

} // namespace arb
