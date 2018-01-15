#pragma once

#include <fstream>
#include <sstream>
#include <math.h>
#include <regex>
#include <tuple>
#include <vector>
#include <utility>


#include <common_types.hpp>
#include <util/debug.hpp>
#include <json/json.hpp>
#include "connection_generator.hpp"

namespace con_gen_util {
    // Simple exception class that can be instantiated with a string as what msg
    class con_gen_error : public std::exception
    {
    public:
        con_gen_error(std::string what) : what_(what) {}

    private:
        virtual const char* what() const throw() {
            return what_.c_str();
        }

    private:
        std::string what_;
    };

    arb::cell_kind cell_kind_from_string(std::string str);


    // Simple population parser. Expects a json like file of populations.
    // On error parsing will stop and what is parsed until then is returned
    // throws a con_gen_errors when there are problems opening or parsing
    // {"population_1":
    //      {
    //          "x_dim":10,
    //           "y_dim" : 10,
    //          "periodic_border" : true,
    //          "cell_type" : "cable1d_neuron"
    //      }
    // }
    std::vector<arb_con_gen::population> parse_populations_from_path(std::string path) {
        std::vector<arb_con_gen::population> populations;

        std::ifstream fid(path);
        if (fid) {
            try {
                nlohmann::json fopts;
                fid >> fopts;

                // Loop over the 'json' population entries
                for (nlohmann::json::iterator it = fopts.begin(); it != fopts.end(); ++it) {
                    try {
                        std::string name = it.key();
                        arb::cell_size_type x_dim = it.value()["x_dim"];
                        arb::cell_size_type y_dim = it.value()["y_dim"];
                        bool periodic = it.value()["periodic_border"];
                        arb::cell_kind kind = cell_kind_from_string(it.value()["cell_type"]);

                        nlohmann::json cell_opts = it.value()["cell_opts"];

                        populations.push_back({name, x_dim, y_dim, periodic, kind, cell_opts});
                    }
                    catch (std::exception& e) {
                        std::cerr << "Could not parse:\n" << it.value() << "\n";
                        throw con_gen_error(
                            "Could not parse entry in: " + path + ": " + e.what());
                    }
                }
            }
            catch (std::exception& e) {
                throw con_gen_error(
                    "unable to parse parameters in " + path + ": " + e.what());
            }
        }
        else {
            throw con_gen_error("Unable to open file" + path);
        }

        return populations;
    }


    // Simple projection parser. Expects a csv like file of population.
    // On error parsing will stop and what is parsed until then is returned
    // throws a con_gen_error when file cannot be opened
    // The lines are parsed separated by a single comma:
    // [idx_pre_polation, idx_post_polation, n_synapses, sd_distance_dep_probabilty
    //  mean_weight, sd_weight, min_delay, delay_per_sd_distance]
    // with types:
    // unsigned, unsigned, unsigned, float, float, float, float, float
    std::vector<arb_con_gen::projection> parse_projections_from_path(std::string path) {

        //unsigned pre_population_id;
        //unsigned post_population_id;
        //unsigned count;
        //float    sd;
        //float    mean_weight;
        //float    sd_weight;
        //float    min_delay;
        //float    delay_per_sd;
        //char comma;

        std::vector<arb_con_gen::projection> projection;

        std::ifstream fid(path);
        if (fid) {
            try {
                nlohmann::json fopts;
                fid >> fopts;

                // Loop over the 'json' population entries
                for (nlohmann::json::iterator it = fopts.begin(); it != fopts.end(); ++it) {
                    try {
                        std::string name = it.key();

                        std::string population_pre = it.value()["population_pre"];
                        std::string population_post = it.value()["population_post"];;
                        unsigned count = it.value()["count"];
                        float    std_2d_kernel = it.value()["std_2d_kernel"];
                        float    weight_mean = it.value()["weight_mean"];
                        float    weight_std = it.value()["weight_std"];
                        float    delay_min = it.value()["delay_min"];
                        float    delay_per_std = it.value()["delay_per_std"];

                        projection.push_back({ population_pre , population_post,
                        {count, std_2d_kernel, weight_mean, weight_std, delay_min, delay_per_std } });

                    }
                    catch (std::exception& e) {
                        std::cerr << "Could not parse:\n" << it.value() << "\n";
                        throw con_gen_error(
                            "Could not parse entry in: " + path + ": " + e.what());
                    }
                }
            }
            catch (std::exception& e) {
                throw con_gen_error(
                    "unable to parse parameters in " + path + ": " + e.what());
            }
        }
        else {
            throw con_gen_error("Unable to open file" + path);
        }

        return projection;
    }

    // Simple gid parser. Expects a comma separated list of individual gid or two gids
    // on a line with the begin and start of a range
    // On error parsing will stop and what is parsed until then is returned
    // throws a con_gen_error when file cannot be opened
    // IF a line starts with a comma, it is parsed as a comma separated list of gids
    // finished with a '<' character
    // If a line starts with a - (or any other character). It is parsed as two
    // comma separated gids and assumed to be a range
    // Types of parsed numbers is unsigned
    std::vector<arb::cell_gid_type> parse_gids_from_path(std::string path) {
        std::ifstream infile(path);

        arb::cell_gid_type gid;
        arb::cell_gid_type gid_untill;
        char char_parsed;
        char comma;

        std::vector<arb::cell_gid_type> gids;

        if (infile) {
            std::string line;
            while (std::getline(infile, line)) {
                std::istringstream iss(line);

                iss >> char_parsed;

                // Parse the first character if it is a , we expect a list of
                // comma separated gids
                if (char_parsed == ',') {
                    while (iss >> gid >> char_parsed) {
                        gids.push_back(gid);
                    }
                }
                // assume it is a pair of gids defining a start and end point
                else {
                    iss >> gid >> comma >> gid_untill;
                    // No error checking!!
                    for (; gid < gid_untill; ++gid) {
                        gids.push_back(gid);
                    }
                }

            }
        }
        else {
            throw con_gen_error("Could not open supplied gids config");
        }

        return gids;
    }

    // Default populations:
    // Two 2d sheets of 100 by 100 neurons with periodic bordes
    std::vector<arb_con_gen::population> default_populations() {
        std::vector<arb_con_gen::population> default_populations;

        default_populations.push_back({"population_1", 10, 10, true, arb::cell_kind::cable1d_neuron });
        default_populations.push_back({"population_2", 10, 10, true, arb::cell_kind::cable1d_neuron });

        return default_populations;
    }

    // Default Gids:
    // Selected such that they are on the border of the sheet. This illustrates
    // The periodic border optimally. Most gids are paired such they the
    // presynaptic neurons arbors overlap
    // 15070, 5030 are shiften in relationship to each other.
    std::vector<arb::cell_gid_type> default_gids() {
        return { 10320, 12003, 17997, 19580,
            15070, 5030,  // These two are shifted !!
            320, 2003, 7997, 9580, 5500 };
    }

    // Default connectome
    // #1:
    // 0 > 1. count 100, ds 0.02 | weight -mean 20.0 -sd 1.0 | delay 1.0 -sd 1.0
    // #2
    // 1 > 0. count 1000, ds 0.1 | weight -mean 2.0 -sd 1.0 | delay 1.0 -sd 1.0
    std::vector<arb_con_gen::projection> default_connectome() {
        std::vector<arb_con_gen::projection>  connectome;
        connectome.push_back({ "population_1","population_2",{  8, 0.02, 2.0, 1.0, 1.0, 1.0 } });
        connectome.push_back({ "population_2","population_1",{  10,0.05, 2.0, 1.0, 1.0, 1.0 } });

        return connectome;
    }

    // Small helper function that converts a string to cell_kind
    arb::cell_kind cell_kind_from_string(std::string str)
    {
        if (str == std::string("cable1d_neuron")) {
            return arb::cell_kind::cable1d_neuron;
        }
        else if (str == std::string("regular_spike_source")) {
            return arb::cell_kind::regular_spike_source;
        }
        else if (str == std::string("data_spike_source")) {
            return arb::cell_kind::data_spike_source;
        }
        else if (str == std::string("inhomogeneous_poisson_spike_source")) {
            return arb::cell_kind::inhomogeneous_poisson_spike_source;
        }
        else
        {
            throw con_gen_error("Unkown cell kind representation encountered: " + str);
        }


    }

}