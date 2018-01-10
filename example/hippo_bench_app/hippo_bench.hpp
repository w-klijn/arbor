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

#include "../con_gen/connection_generator.hpp"

namespace hippo_bench {
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


    // Simple population parser. Expects a csv like file of population.
    // On error parsing will stop and what is parsed until then is returned
    // throws a con_gen_error when file cannot be opened
    // The lines are parsed separated by a single comma:
    // [x_dim, y_dim, periodic]
    // with types:
    // unsigned, unsigned, 0 -or- 1
    std::vector<arb_con_gen::population> parse_populations_from_path(std::string path) {
        std::ifstream infile(path);
        arb::cell_gid_type x_dim;
        arb::cell_gid_type y_dim;
        bool periodic;
        char comma;

        std::vector<arb_con_gen::population> populations;

        // Regex parsing start
        //std::regex e;
        //std::smatch match;
        //e = "^\\s*(\\d+)\\s*,\\s*(\\d+)\\s*,\\s*((?i)(true|false))\\s*";  // usigned , usigned , true
        //std::cout << "debug" << std::endl;
        //if (infile) {
        //    std::string line;

        //    while (std::getline(infile, line)) {
        //        std::istringstream iss(line);

        //        std::cout << "debug" << std::endl;
        //        if (std::regex_search(line, match, e))
        //        {
        //            std::cout << "match: " << match[0] << '\n';
        //            std::cout << "match: " << match[1] << '\n';
        //            std::cout << "match: " << match[2] << '\n';

        //        }


        //        populations.push_back({ x_dim, y_dim, periodic });

        //    }
        //}


        if (infile) {
            std::string line;
            while (std::getline(infile, line)){
                std::istringstream iss(line);
                if (!(iss >> x_dim >> comma >> y_dim >> comma >> periodic)) {
                    break; }
                populations.push_back({ x_dim, y_dim, periodic });

            }
        }
        else {
            throw con_gen_error("Could not open supplied population config");
        }
        std::cout << "Debug 2" << std::endl;
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
        std::ifstream infile(path);

        unsigned pre_population_id;
        unsigned post_population_id;
        unsigned count;
        float    sd;
        float    mean_weight;
        float    sd_weight;
        float    min_delay;
        float    delay_per_sd;
        char comma;

        std::vector<arb_con_gen::projection> projection;

        if (infile) {
            std::string line;
            while (std::getline(infile, line)) {
                std::istringstream iss(line);

                if (!(iss >> pre_population_id >> comma >> post_population_id >> comma >>
                count >> comma >> sd >> comma >> mean_weight >> comma >>
                sd_weight >> comma >> min_delay >> comma >> delay_per_sd)) {
                    break;
                }
                projection.push_back({ pre_population_id,post_population_id,{
                    sd, count, mean_weight, sd_weight,min_delay, delay_per_sd}});
            }
        }
        else {
            throw con_gen_error("Could not open supplied projection config");
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

        default_populations.push_back({ 100, 100, true });
        default_populations.push_back({ 100, 100, true });

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
        connectome.push_back({ 0,1,{ 0.02, 400,  2.0, 1.0, 1.0, 1.0 } });
        connectome.push_back({ 1,0,{ 0.05, 1000, 2.0, 1.0, 1.0, 1.0 } });

        return connectome;
    }

}