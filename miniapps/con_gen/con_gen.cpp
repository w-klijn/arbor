#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <utility>

#include <tinyopt.hpp>
#include <util/optional.hpp>
#include <common_types.hpp>

#include "connection_generator.hpp"


namespace to = arb::to;
using arb::util::optional;

const char* usage_str =
"[OPTION]...\n"
"\n"
"  -n, --count=int        (10000)  Number of individual Poisson cell to run.\n"
"  -b, --begin=float      (0.0)    Start time (in ms) when to start generating\n"
"  -e, --end=float        (1000.0) End time (in ms) when to end generation\n"
"  -s, --sample=float     (0.1)    Internal sample rate for Poisson process\n"
"  -i, --interpolate=0/1  (1)   Interpolate between the supplied time-rate pairs\n"
"\n"
"  --pairs=path           (option)     Path to file with 'float, float' time rate pairs\n"
"  --output=path          (./spikes.gdf) Export produced spikes to this path\n"

"  -h, --help             Emit this message and exit.\n"
"\n"
"Create a group of Inhomogeneous Poisson Spike Sources\n"
"Run the cells and output the produced spikes to file\n"
"\n"
"Output of the default settings can be used in parse_and_plot.py to generate\n"
"the example plot in the how-to at 'http://arbor.readthedocs.io/en/latest/'\n"
"\n"
" Default time varying inhomogeneous spike rate:\n"
" hz.|                     \n"
" 240|     _-_             \n"
"    |    -   -  -         \n"
"    |   -     -- -        \n"
" 0  |__-__________-__     \n"
"      100        900   ms \n"
;

int main(int argc, char** argv) {
    // options
    //unsigned n_cells = 10000;
    //arb::time_type begin = 0.0;
    //arb::time_type end = 1000.0;
    //arb::time_type sample_delta = 0.1;
    //bool interpolate = true;

    try {
        auto arg = argv + 1;
        while (*arg) {
            //if (auto o = to::parse_opt<int>(arg, 'n', "count")) {
            //    n_cells = *o;
            //}
            //else if (auto o = to::parse_opt<arb::time_type>(arg, 'b', "begin")) {
            //    begin = *o;
            //}
            //if
            //else {
                throw to::parse_opt_error(*arg, "unrecognized option");
            //}
        }
    }
    catch (to::parse_opt_error& e) {
        std::cerr << argv[0] << ": " << e.what() << "\n";
        std::cerr << "Try '" << argv[0] << " --help' for more information.\n";
        std::exit(2);
    }
    catch (std::exception& e) {
        std::cerr << "caught exception: " << e.what() << "\n";
        std::exit(1);
    }

    // Create two population of 100 by 100 cells
    std::vector<arb::population> populations;
    populations.push_back({ 100, 100, true  });
    populations.push_back({ 100, 100, true });

    // Create a projection from index 0 to index 1
    std::vector<std::tuple<unsigned, unsigned, arb::projection_pars>>  connectome;
    arb::projection_pars pars(0.02, 10, 2.0, 1.0, 1.0, 1.0);
    std::tuple<unsigned, unsigned, arb::projection_pars> proj(0, 1, pars);
    connectome.push_back(proj);
    arb::projection_pars pars2(0.1, 10, 2.0, 1.0, 1.0, 1.0);
    std::tuple<unsigned, unsigned, arb::projection_pars> proj2(0, 1, pars2);
    connectome.push_back(proj2);



    arb::connection_generator gen(populations, connectome);
    //std::vector<arb::cell_gid_type> gids = gen.pre_synaptic_cells(15050);

    auto synapses = gen.synapses_on(10099);

    std::ofstream outfile("gids.dat");
    if (outfile) {
        //for (auto gid : gids) {
        //    outfile << gid << "\n";
        //}

        for (auto synapse : synapses) {
            outfile << synapse.gid << "," << synapse.weight << "," << synapse.delay << "\n";
        }

    }


    return 0;
}
