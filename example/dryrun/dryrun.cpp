/*
 * A miniapp that demonstrates how to use dry_run mode
 *
 */

#include <fstream>
#include <iomanip>
#include <iostream>

#include <nlohmann/json.hpp>

#include <arbor/assert_macro.hpp>
#include <arbor/common_types.hpp>
#include <arbor/context.hpp>
#include <arbor/load_balance.hpp>
#include <arbor/cable_cell.hpp>
#include <arbor/profile/meter_manager.hpp>
#include <arbor/profile/profiler.hpp>
#include <arbor/simple_sampler.hpp>
#include <arbor/simulation.hpp>
#include <arbor/symmetric_recipe.hpp>
#include <arbor/recipe.hpp>
#include <arbor/version.hpp>

#include <sup/ioutil.hpp>
#include <sup/json_meter.hpp>
#include <sup/json_params.hpp>

#ifdef ARB_MPI_ENABLED
#include <mpi.h>
#include <arborenv/with_mpi.hpp>
#endif

// Parameters used to generate the random cell morphologies.
struct cell_parameters {
    //  Maximum number of levels in the cell (not including the soma)
    unsigned max_depth = 5;

    // The following parameters are described as ranges.
    // The first value is at the soma, and the last value is used on the last level.
    // Values at levels in between are found by linear interpolation.
    std::array<double,2> branch_probs = {1.0, 0.5}; //  Probability of a branch occuring.
    std::array<unsigned,2> compartments = {20, 2};  //  Compartment count on a branch.
    std::array<double,2> lengths = {200, 20};       //  Length of branch in μm.
};

struct run_params {
    std::string name = "default";
    bool dry_run = false;
    unsigned num_cells_per_rank = 10;
    unsigned num_ranks = 1;
    double min_delay = 10;
    double duration = 100;
    cell_parameters cell;
};

run_params read_options(int argc, char** argv);

using arb::cell_gid_type;
using arb::cell_lid_type;
using arb::cell_size_type;
using arb::cell_member_type;
using arb::cell_kind;
using arb::time_type;
using arb::cell_probe_address;

// Generate a cell.
arb::cable_cell branch_cell(arb::cell_gid_type gid, const cell_parameters& params);

class tile_desc: public arb::tile {
public:
    tile_desc(unsigned num_cells, unsigned num_tiles, cell_parameters params, unsigned min_delay):
            num_cells_(num_cells),
            num_tiles_(num_tiles),
            cell_params_(params),
            min_delay_(min_delay)
    {}

    cell_size_type num_cells() const override {
        return num_cells_;
    }

    cell_size_type num_tiles() const override {
        return num_tiles_;
    }

    arb::util::unique_any get_cell_description(cell_gid_type gid) const override {
        return branch_cell(gid, cell_params_);
    }

    cell_kind get_cell_kind(cell_gid_type gid) const override {
        return cell_kind::cable;
    }

    arb::util::any get_global_properties(arb::cell_kind) const override {
        arb::cable_cell_global_properties gprop;
        gprop.default_parameters = arb::neuron_parameter_defaults;
        return gprop;
    }

    // Each cell has one spike detector (at the soma).
    cell_size_type num_sources(cell_gid_type gid) const override {
        return 1;
    }

    // The cell has one target synapse
    cell_size_type num_targets(cell_gid_type gid) const override {
        return 1;
    }

    // Each cell has one incoming connection, from any cell in the network spanning all ranks:
    // src gid in {0, ..., num_cells_*num_tiles_ - 1}.
    std::vector<arb::cell_connection> connections_on(cell_gid_type gid) const override {
        std::uniform_int_distribution<cell_gid_type> source_distribution(0, num_cells_*num_tiles_ - 2);

        auto src_gen = std::mt19937(gid);
        auto src = source_distribution(src_gen);
        if (src>=gid) ++src;

        return {arb::cell_connection({src, 0}, {gid, 0}, event_weight_, min_delay_)};
    }

    // Return an event generator on every 20th gid. This function needs to generate events
    // for ALL cells on ALL ranks. This is because the symmetric recipe can not easily
    // translate the src gid of an event generator
    std::vector<arb::event_generator> event_generators(cell_gid_type gid) const override {
        std::vector<arb::event_generator> gens;
        if (gid%20 == 0) {
            gens.push_back(arb::explicit_generator(arb::pse_vector{{{gid, 0}, 0.1, 1.0}}));
        }
        return gens;
    }

    // There is one probe (for measuring voltage at the soma) on the cell.
    cell_size_type num_probes(cell_gid_type gid)  const override {
        return 1;
    }

    arb::probe_info get_probe(cell_member_type id) const override {
        // Get the appropriate kind for measuring voltage.
        cell_probe_address::probe_kind kind = cell_probe_address::membrane_voltage;
        // Measure at the soma.
        arb::segment_location loc(0, 0.0);

        return arb::probe_info{id, kind, cell_probe_address{loc, kind}};
    }

private:
    cell_size_type num_cells_;
    cell_size_type num_tiles_;
    cell_parameters cell_params_;
    double min_delay_;
    float event_weight_ = 0.01;
};

struct cell_stats {
    using size_type = unsigned;
    size_type ncells = 0;
    int nranks = 1;
    size_type nsegs = 0;
    size_type ncomp = 0;

    cell_stats(arb::recipe& r, run_params params) {
#ifdef ARB_MPI_ENABLED
        if(!params.dry_run) {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &nranks);
            ncells = r.num_cells();
            size_type cells_per_rank = ncells/nranks;
            size_type b = rank*cells_per_rank;
            size_type e = (rank+1)*cells_per_rank;
            size_type nsegs_tmp = 0;
            size_type ncomp_tmp = 0;
            for (size_type i=b; i<e; ++i) {
                auto c = arb::util::any_cast<arb::cable_cell>(r.get_cell_description(i));
                nsegs_tmp += c.num_segments();
                ncomp_tmp += c.num_compartments();
            }
            MPI_Allreduce(&nsegs_tmp, &nsegs, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&ncomp_tmp, &ncomp, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
        }
#else
        if(!params.dry_run) {
            nranks = 1;
            ncells = r.num_cells();
            for (size_type i = 0; i < ncells; ++i) {
                auto c = arb::util::any_cast<arb::cable_cell>(r.get_cell_description(i));
                nsegs += c.num_segments();
                ncomp += c.num_compartments();
            }
        }
#endif
        else {
            nranks = params.num_ranks;
            ncells = r.num_cells(); //total number of cells across all ranks

            for (size_type i = 0; i < params.num_cells_per_rank; ++i) {
                auto c = arb::util::any_cast<arb::cable_cell>(r.get_cell_description(i));
                nsegs += c.num_segments();
                ncomp += c.num_compartments();
            }

            nsegs *= params.num_ranks;
            ncomp *= params.num_ranks;
        }
    }

    friend std::ostream& operator<<(std::ostream& o, const cell_stats& s) {
        return o << "cell stats: "
                 << s.nranks << " ranks; "
                 << s.ncells << " cells; "
                 << s.nsegs << " segments; "
                 << s.ncomp << " compartments.";
    }
};


int main(int argc, char** argv) {
    try {
#ifdef ARB_MPI_ENABLED
        arbenv::with_mpi guard(argc, argv, false);
#endif
        bool root = true;
        auto params = read_options(argc, argv);

        auto resources = arb::proc_allocation();
        auto ctx = arb::make_context(resources);

        if (params.dry_run) {
            ctx = arb::make_context(resources, arb::dry_run_info(params.num_ranks, params.num_cells_per_rank));
        }
#ifdef ARB_MPI_ENABLED
        else {
            ctx = arb::make_context(resources, MPI_COMM_WORLD);
            {
                int rank;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                root = rank==0;
            }
        }
#endif
        arb_assert(arb::num_ranks(ctx)==params.num_ranks);


#ifdef ARB_PROFILE_ENABLED
        arb::profile::profiler_initialize(ctx);
#endif
        std::cout << sup::mask_stream(root);

        // Print a banner with information about hardware configuration
        std::cout << "gpu:      " << (has_gpu(ctx)? "yes": "no") << "\n";
        std::cout << "threads:  " << num_threads(ctx) << "\n";
        std::cout << "mpi:      " << (has_mpi(ctx)? "yes": "no") << "\n";
        std::cout << "ranks:    " << num_ranks(ctx) << "\n" << std::endl;

        std::cout << "run mode: " << distribution_type(ctx) << "\n";

        arb::profile::meter_manager meters;
        meters.start(ctx);

        // Create an instance of our tile and use it to make a symmetric_recipe.
        auto tile = std::make_unique<tile_desc>(params.num_cells_per_rank,
                params.num_ranks, params.cell, params.min_delay);
        arb::symmetric_recipe recipe(std::move(tile));

        cell_stats stats(recipe, params);
        std::cout << stats << "\n";

        auto decomp = arb::partition_load_balance(recipe, ctx);

        // Construct the model.
        arb::simulation sim(recipe, decomp, ctx);

        // Set up recording of spikes to a vector on the root process.
        std::vector<arb::spike> recorded_spikes;
        if (root) {
            sim.set_global_spike_callback(
                    [&recorded_spikes](const std::vector<arb::spike>& spikes) {
                        recorded_spikes.insert(recorded_spikes.end(), spikes.begin(), spikes.end());
                    });
        }

        meters.checkpoint("model-init", ctx);

        // Run the simulation for 100 ms, with time steps of 0.025 ms.
        sim.run(params.duration, 0.025);

        meters.checkpoint("model-run", ctx);

        auto ns = sim.num_spikes();
        std::cout << "\n" << ns << " spikes generated at rate of "
                  << params.duration/ns << " ms between spikes\n\n";

        // Write spikes to file
        if (root) {
            std::ofstream fid("spikes.gdf");
            if (!fid.good()) {
                std::cerr << "Warning: unable to open file spikes.gdf for spike output\n";
            }
            else {
                char linebuf[45];
                for (auto spike: recorded_spikes) {
                    auto n = std::snprintf(
                            linebuf, sizeof(linebuf), "%u %.4f\n",
                            unsigned{spike.source.gid}, float(spike.time));
                    fid.write(linebuf, n);
                }
            }
        }

        auto profile = arb::profile::profiler_summary();
        std::cout << profile << "\n";

        auto report = arb::profile::make_meter_report(meters, ctx);
        std::cout << report;
    }
    catch (std::exception& e) {
        std::cerr << "exception caught in ring miniapp:\n" << e.what() << "\n";
        return 1;
    }

    return 0;
}

// Helper used to interpolate in branch_cell.
template <typename T>
double interp(const std::array<T,2>& r, unsigned i, unsigned n) {
    double p = i * 1./(n-1);
    double r0 = r[0];
    double r1 = r[1];
    return r[0] + p*(r1-r0);
}

arb::cable_cell branch_cell(arb::cell_gid_type gid, const cell_parameters& params) {
    arb::cable_cell cell;
    cell.default_parameters.axial_resistivity = 100; // [Ω·cm]

    // Add soma.
    auto soma = cell.add_soma(12.6157/2.0); // For area of 500 μm².
    soma->add_mechanism("hh");

    std::vector<std::vector<unsigned>> levels;
    levels.push_back({0});

    // Standard mersenne_twister_engine seeded with gid.
    std::mt19937 gen(gid);
    std::uniform_real_distribution<double> dis(0, 1);

    double dend_radius = 0.5; // Diameter of 1 μm for each cable.

    unsigned nsec = 1;
    for (unsigned i=0; i<params.max_depth; ++i) {
        // Branch prob at this level.
        double bp = interp(params.branch_probs, i, params.max_depth);
        // Length at this level.
        double l = interp(params.lengths, i, params.max_depth);
        // Number of compartments at this level.
        unsigned nc = std::round(interp(params.compartments, i, params.max_depth));

        std::vector<unsigned> sec_ids;
        for (unsigned sec: levels[i]) {
            for (unsigned j=0; j<2; ++j) {
                if (dis(gen)<bp) {
                    sec_ids.push_back(nsec++);
                    auto dend = cell.add_cable(sec, arb::section_kind::dendrite, dend_radius, dend_radius, l);
                    dend->set_compartments(nc);
                    dend->add_mechanism("pas");
                }
            }
        }
        if (sec_ids.empty()) {
            break;
        }
        levels.push_back(sec_ids);
    }

    // Add spike threshold detector at the soma.
    cell.add_detector({0,0}, 10);

    // Add a synapse to the mid point of the first dendrite.
    cell.add_synapse({1, 0.5}, "expsyn");

    return cell;
}

run_params read_options(int argc, char** argv) {
    using sup::param_from_json;

    run_params params;
    if (argc<2) {
        std::cout << "Using default parameters.\n";
        return params;
    }
    if (argc>2) {
        throw std::runtime_error("More than command line one option not permitted.");
    }

    std::string fname = argv[1];
    std::cout << "Loading parameters from file: " << fname << "\n";
    std::ifstream f(fname);

    if (!f.good()) {
        throw std::runtime_error("Unable to open input parameter file: "+fname);
    }

    nlohmann::json json;
    json << f;

    param_from_json(params.name, "name", json);
    param_from_json(params.dry_run, "dry-run", json);
    param_from_json(params.num_cells_per_rank, "num-cells-per-rank", json);
    param_from_json(params.num_ranks, "num-ranks", json);
    param_from_json(params.duration, "duration", json);
    param_from_json(params.min_delay, "min-delay", json);
    param_from_json(params.cell.max_depth, "depth", json);
    param_from_json(params.cell.branch_probs, "branch-probs", json);
    param_from_json(params.cell.compartments, "compartments", json);
    param_from_json(params.cell.lengths, "lengths", json);

    if (!json.empty()) {
        for (auto it=json.begin(); it!=json.end(); ++it) {
            std::cout << "  Warning: unused input parameter: \"" << it.key() << "\"\n";
        }
        std::cout << "\n";
    }

    return params;
}
