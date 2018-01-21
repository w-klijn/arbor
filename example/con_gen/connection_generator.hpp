//#include <domain_decomposition.hpp>
//#include <hardware/node_info.hpp>
//#include <recipe.hpp>

#pragma once

#include <vector>
#include <utility>
#include <tuple>
#include <math.h>
#include <map>
#include <string>


#include <common_types.hpp>
#include <math.hpp>
#include <random>
#include <util/debug.hpp>
#include <json/json.hpp>
#include <recipe.hpp>
#include <util/unique_any.hpp>

namespace arb_con_gen {

// Describes a 2d surface of neurons located on grid locations
// -x_dim   number of neurons on the x-side
// -y_dim   number of neurons on the y-side
// -periodic Do the border loop back to the other side (torus topology)
struct population {
    std::string name;
    arb::cell_size_type x_dim;
    arb::cell_size_type y_dim;
    bool periodic;

    arb::cell_size_type n_cells;
    arb::cell_kind kind;
    nlohmann::json cell_opts;

    // TODO: enum topology_type ( grid, pure random, minimal distance)

    population(std::string name, arb::cell_size_type x_dim, arb::cell_size_type y_dim, bool per,
        arb::cell_kind kind, nlohmann::json cell_opts = {}) :
        name(name), x_dim(x_dim), y_dim(y_dim), periodic(per), n_cells(x_dim *y_dim), kind(kind), cell_opts(cell_opts)
    {

        // Sanity check
        EXPECTS(x_dim > 0);
        EXPECTS(y_dim > 0);
    }


};

// Describes a projection between the neurons between two populations
// - sd      sd of the normal distributed used to sample the pre_synaptic
//           The dimensions of the pre-population is sampled as if it has size 1.0 * 1.0
// - count   Number of samples to take. When sampling from a non periodic population
//           this count can be lower (akin with a sample in-vitro)
//
// - weight_mean  Mean synaptic weight for the created synapse
// - weight_sd    Standard deviation around mean for sampling the weights
//
// - delay_min      Minimal delay of the created synapse
// - delay_per_sd   Delay increase by sd distance between neurons
struct projection_pars {
    arb::cell_size_type count=0;
    float sd = 0.0;

    // parameters for the synapses on this projection
    float weight_mean=0;
    float weight_sd=0;

    float delay_min=0;        // Minimal delay
    float delay_per_sd=0;    // per

    projection_pars( arb::cell_size_type count, float sd, float weight_mean,
        float weight_std, float delay_min, float delay_per_std) :
        count(count), sd(sd),
        weight_mean(weight_mean), weight_sd(weight_std),
        delay_min(delay_min), delay_per_sd(delay_per_std)
    {
        // Sanity checks
        EXPECTS(sd > 0.0);
        EXPECTS(count > 0);
        EXPECTS(delay_min > 0.9999); // TODO: This a neuroscientific 'fact' not needed for valid functioning
        EXPECTS(delay_per_std > 0.0);
    }
};

// Helper struct to collect some parameters together for creating projections
// -pre_name    The index in the population list that is pre synaptic for this
//             projection
// -post_name   The index in the population list that is post synaptic for this
//             projection
// -pars       Parameters used to generate the synapses for this connection
struct projection {
    std::string pre_name;
    std::string post_name;
    projection_pars pars;

    projection(std::string pre_population, std::string post_population, projection_pars pars) :
        pre_name(pre_population), post_name(post_population), pars(pars)
    {}
};

// Return type for connection generation
// A set of pre-synaptic cell gid,
// weight and delay
struct synaps_pars {
    arb::cell_gid_type gid;
    float weight;
    float delay;
    synaps_pars(arb::cell_gid_type gid, float weight, float delay) :
        gid(gid), weight(weight), delay(delay)
    {}
};


// Data structure collecting the parameters needed for Poisson event generator
struct poisson_event_pars {
    double rate;
    double weight;
    double start;

    poisson_event_pars(double rate, double weight, double start) :
        rate(rate), weight(weight), start(start)
    {}
};


// Data structure for collecting all parameters for a cell
struct cell_pars {
    unsigned compartments_per_segment;
    std::string synapse_type;
    std::string dendrite_mechanism;
    double dendrite_rL;
    std::string soma_mechanism;
    unsigned synapses_per_cell;
    double spike_threshold;


    cell_pars(unsigned compartments_per_segment, std::string synapse_type,
        std::string dendrite_mechanism, double dendrite_rL, std::string soma_mechanism,
        unsigned synapses_per_cell, double spike_threshold) :
        compartments_per_segment(compartments_per_segment),
        synapse_type(synapse_type),
        dendrite_mechanism(dendrite_mechanism),
        dendrite_rL(dendrite_rL),
        soma_mechanism(soma_mechanism),
        synapses_per_cell(synapses_per_cell),
        spike_threshold(spike_threshold)
    {}

    // Needed to suppress warnings, but should never be called
    cell_pars() {
        EXPECTS(false);
    }
};

class connection_generator {

public:
    // Connection generator.
    // Expects a vector of populations descriptions and vector of projections
    // between these
    // TODO: This is a first implementation: sub populations are NOT implemented
    connection_generator(const std::vector<population> & populations,
        std::vector<projection> connectome) :
        connectome_(std::move(connectome))
    {
        arb::cell_gid_type gid_idx = 0;

        // Create the local populations with start index set
        for (auto pop : populations) {

            populations_.insert({ pop.name, {pop, gid_idx } });
            gid_idx += pop.n_cells;
        }

        n_cells_ = gid_idx;
    }

    connection_generator() {
        n_cells_ = 0;
    }


    // Get the total count of cells on this connection generator
    // TODO: Refactor to connection generator and neuron generator
    arb::cell_size_type num_cells() const {
        return n_cells_;
    }

    arb::cell_kind get_cell_kind(arb::cell_gid_type gid) const {

        EXPECTS(gid < n_cells_);
        arb::cell_kind kind;
        for (const auto& pop: populations_) {
            if (gid >= pop.second.start_index && gid < pop.second.end_index)
                kind = pop.second.kind;
        }

        return kind;
    }


    // Returns a struct with the cell parameters
    arb::util::unique_any get_cell_opts(arb::cell_gid_type gid) const {
        EXPECTS(gid < n_cells_);
        for (const auto& pop : populations_) {
            if (gid >= pop.second.start_index && gid < pop.second.end_index) {
                auto const& cell_pars_json = pop.second.cell_opts;

                cell_pars pars(
                    cell_pars_json["compartments_per_segment"],
                    cell_pars_json["synapse_type"],
                    cell_pars_json["dendrite_mechanism"],
                    cell_pars_json["dendrite_rL"],
                    cell_pars_json["soma_mechanism"],
                    cell_pars_json["synapses_per_cell"],
                    cell_pars_json["spike_threshold"]);


                return arb::util::unique_any(std::move(pars));
            }
        }

        EXPECTS(false);
        return {};

    }


    // Return the poisson_generators on this cell
    std::vector<arb_con_gen::poisson_event_pars> const get_cell_poisson_generators(arb::cell_gid_type gid) const {
        EXPECTS(gid < n_cells_);

        std::vector<arb_con_gen::poisson_event_pars> generator_pars;

        for (const auto& pop : populations_) {
            // Find the matching population
            if (gid >= pop.second.start_index && gid < pop.second.end_index) {
                // Json description
                auto const& cell_pars_json = pop.second.cell_opts;

                // If no generators declared
                if (cell_pars_json.find("poisson_generators") == cell_pars_json.end()) {
                     return generator_pars;
                }

                //Loop over the generators
                for (nlohmann::json::const_iterator it = cell_pars_json["poisson_generators"].begin();
                    it != cell_pars_json["poisson_generators"].end(); ++it) {
                    generator_pars.push_back({
                        it.value()["rate"], it.value()["weight"], it.value()["start"]
                    });
                }

                return generator_pars;
            }
        }

        EXPECTS(false);
        return {};
    }


    // Returns the number of synapses on this cell
    arb::cell_size_type num_synapses_on(arb::cell_gid_type gid) const {

        std::mt19937 gen;
        gen.seed(gid);
        unsigned synapse_count = 0;
        for (auto project : connectome_) {
            // Sanity check that the populations exist
            EXPECTS(populations_.count(project.pre_name));
            EXPECTS(populations_.count(project.post_name));

            // Shorthand for the pre and post populations
            auto pre_pop = populations_.at(project.pre_name);
            auto post_pop = populations_.at(project.post_name);
            auto pro_pars = project.pars;

            // Check if this gid receives connections via this projection
            // TODO: Replace with the fancy in range function we have somewhere in the utils
            if (gid < post_pop.start_index || gid >= (post_pop.end_index)) {
                continue;
            }


            // Convert the gid to a index location on a grid
            point post_location = {
                float((gid - post_pop.start_index) % post_pop.x_dim) / post_pop.x_dim,
                float((gid - post_pop.start_index) / post_pop.y_dim) / post_pop.y_dim };

            // If we have non square sides we need to correct the stdev to get
            // circular projections! We do all c
            float sd_x = pro_pars.sd;
            float sd_y = pro_pars.sd;
            if (post_pop.x_dim != post_pop.y_dim) {
                if (post_pop.x_dim < post_pop.y_dim) {
                    float ratio = float(post_pop.y_dim) / post_pop.x_dim;
                    sd_x *= ratio;
                }
                else {
                    float ratio = float(post_pop.x_dim) / post_pop.y_dim;
                    sd_y *= ratio;
                }
            }

            // The distribution for these locations
            std::normal_distribution<float> distr_x(post_location.x, sd_x);
            std::normal_distribution<float> distr_y(post_location.y, sd_y);

            // Start generating the connections
            for (arb::cell_gid_type idx = 0; idx < pro_pars.count; ++idx) {
                // draw the locations
                float x_source = distr_x(gen);
                float y_source = distr_y(gen);
                if (post_pop.periodic) {
                    // Todo: add non-periodic borders
                    // normalize: move all values between [0.0, 1.0)
                    x_source -= int(std::floor(x_source));
                    y_source -= int(std::floor(y_source));
                }
                else {
                    // If we have non periodic borders this connection is not
                    // created (akin to a in vitro slice) if outside of [0, 1.0)
                    if (x_source < 0 || x_source >= 1.0 || y_source < 0 || y_source >= 1.0) {
                        continue;
                    }
                }
                synapse_count++;
            }
        }

        return synapse_count;
    }

    // Returns a vector of all synaptic parameters sets for this gid
    std::vector<arb::cell_connection> synapses_on(arb::cell_gid_type gid) const {
        std::mt19937 gen;
        gen.seed(gid);

        std::vector<arb::cell_connection> connections;
        unsigned idx_projection = 0;
        for (auto project : connectome_) {
            // Sanity check that the populations exist
            EXPECTS(populations_.count(project.pre_name));
            EXPECTS(populations_.count(project.post_name));

            // Shorthand for the pre and post populations
            auto pre_pop = populations_.at(project.pre_name);
            auto post_pop = populations_.at(project.post_name);
            auto pro_pars = project.pars;

            // Check if this gid receives connections via this projection
            // TODO: Replace with the fancy in range function we have somewhere in the utils
            if (gid < post_pop.start_index || gid >= (post_pop.end_index)) {
                continue;
            }

            // Distribution to draw the weights
            std::normal_distribution<float> weight_distr(pro_pars.weight_mean, pro_pars.weight_sd);

            // Weight sign
            float weight_sign = arb::math::signum(pro_pars.weight_mean);

            // Convert the gid to a index location on a grid
            point post_location = {
                float((gid - post_pop.start_index) % post_pop.x_dim) / post_pop.x_dim,
                float((gid - post_pop.start_index) / post_pop.y_dim) / post_pop.y_dim};

            // If we have non square sides we need to correct the stdev to get
            // circular projections! We do all c
            float sd_x = pro_pars.sd;
            float sd_y = pro_pars.sd;
            if (post_pop.x_dim != post_pop.y_dim) {
                if (post_pop.x_dim < post_pop.y_dim) {
                    float ratio = float(post_pop.y_dim) / post_pop.x_dim;
                    sd_x *= ratio;
                }
                else {
                    float ratio = float(post_pop.x_dim) / post_pop.y_dim;
                    sd_y *= ratio;
                }
            }

            // The distribution for these locations
            std::normal_distribution<float> distr_x(post_location.x, sd_x);
            std::normal_distribution<float> distr_y(post_location.y, sd_y);
            float mean_sd = (sd_x + sd_y) / 2.0;

            // Start generating the connections
            for (arb::cell_gid_type idx = 0; idx < pro_pars.count; ++idx) {
                // draw the locations
                float x_source = distr_x(gen);
                float y_source = distr_y(gen);

                // We need the distance between the points to calculate the delay.
                // This is calculated easy before periodization of the locations
                // TODOW: With the two sd different this
                float weighted_distance = std::sqrt(
                    std::pow(post_location.x - x_source, 2) +
                    std::pow(post_location.y - y_source, 2)) / mean_sd;

                if (post_pop.periodic) {
                    // Todo: add non-periodic borders
                    // normalize: move all values between [0.0, 1.0)
                    x_source -= int(std::floor(x_source));
                    y_source -= int(std::floor(y_source));
                }
                else {
                    // If we have non periodic borders this connection is not
                    // created (akin to a in vitro slice) if outside of [0, 1.0)
                    if (x_source < 0 || x_source >= 1.0 || y_source < 0 || y_source >= 1.0) {
                        continue;
                    }
                }

                arb::cell_gid_type gid_pre = arb::cell_gid_type(y_source * pre_pop.y_dim) * pre_pop.x_dim +
                    arb::cell_gid_type(x_source * pre_pop.x_dim);
                // absolute gid
                gid_pre += pre_pop.start_index;

                // TODO: If we have randomly distributed cell, use a quadtree to find the gid

                float delay = weighted_distance * pro_pars.delay_per_sd + pro_pars.delay_min;
                float weight = weight_distr(gen);

                // Flip the sign of the weight depending if we are incorrect
                weight = (weight_sign * weight) < 0?  -weight: weight;

                connections.push_back({ { gid, 0 }, { gid_pre, 0 }, weight, delay });
            }
        }

        return connections;
    }

private:

    struct point {
        float x;
        float y;
    };


    struct population_instantiated : public population {
        arb::cell_gid_type start_index;
        arb::cell_gid_type end_index;

        population_instantiated(population pop, arb::cell_gid_type start_index) :
            population(pop), start_index(start_index), end_index(start_index + pop.n_cells )
        {}
    };

    std::map<std::string, population_instantiated> populations_;
    std::vector<projection> connectome_;

    // Number of cells in this connection class
    arb::cell_size_type n_cells_;

};

}