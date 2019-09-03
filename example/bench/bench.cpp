/*
 * Miniapp that uses the artificial benchmark cell type to test
 * the simulator infrastructure.
 */
#include <fstream>
#include <iomanip>
#include <iostream>

#include <nlohmann/json.hpp>

#include <arbor/profile/meter_manager.hpp>
#include <arbor/benchmark_cell.hpp>
#include <arbor/context.hpp>
#include <arbor/domain_decomposition.hpp>
#include <arbor/load_balance.hpp>
#include <arbor/profile/profiler.hpp>
#include <arbor/recipe.hpp>
#include <arbor/simulation.hpp>
#include <arbor/version.hpp>


#include <arborenv/concurrency.hpp>
#include <arborenv/gpu_env.hpp>

#include <sup/ioutil.hpp>
#include <sup/json_meter.hpp>
#include <sup/json_params.hpp>

#ifdef ARB_MPI_ENABLED
#include <mpi.h>
#include <arborenv/with_mpi.hpp>
#endif

struct bench_params {
    struct cell_params {
        double spike_freq_hz = 10;   // Frequency in hz that cell will generate (poisson) spikes.
        double realtime_ratio = 0.1; // Integration speed relative to real time, e.g. 10 implies
                                     // that a cell is integrated 10 times slower than real time.
    };
    struct network_params {
        unsigned fan_in = 5000;      // Number of incoming connections on each cell.
        double min_delay = 10;       // Used as the delay on all connections.
    };
    std::string name = "default";    // Name of the model.
    unsigned num_cells = 1000;       // Number of cells in model.
    arb::time_type duration = 100;          // Simulation duration in ms.

    cell_params cell;                // Cell parameters for all cells in model.
    network_params network;          // Description of the network.

    // Expected simulation performance properties based on model parameters.

    // Time to finish simulation if only cell overheads are counted.
    double expected_advance_time() const {
        return cell.realtime_ratio * duration*1e-3 * num_cells;
    }
    // Total expected number of spikes generated by simulation.
    unsigned expected_spikes() const {
        return num_cells * duration*1e-3 * cell.spike_freq_hz;
    }
    // Expected number of spikes generated per min_delay/2 interval.
    unsigned expected_spikes_per_interval() const {
        return num_cells * network.min_delay*1e-3/2 * cell.spike_freq_hz;
    }
    // Expected number of post-synaptic events delivered over simulation.
    unsigned expected_events() const {
        return expected_spikes() * network.fan_in;
    }
    // Expected number of post-synaptic events delivered per min_delay/2 interval.
    unsigned expected_events_per_interval() const {
        return expected_spikes_per_interval() * network.fan_in;
    }
};

bench_params read_options(int argc, char** argv);
std::ostream& operator<<(std::ostream& o, const bench_params& p);

class bench_recipe: public arb::recipe {
    bench_params params_;

public:
    bench_recipe(bench_params p): params_(std::move(p)) {}

    arb::cell_size_type num_cells() const override {
        return params_.num_cells;
    }

    arb::util::unique_any get_cell_description(arb::cell_gid_type gid) const override {
        std::mt19937_64 rng(gid);
        arb::benchmark_cell cell;
        cell.realtime_ratio = params_.cell.realtime_ratio;

        // The time_sequence of the cell produces the series of time points at
        // which it will spike. We use a poisson_schedule with a random sequence
        // seeded with the gid. In this way, a cell's random stream depends only
        // on its gid, and will hence give reproducable results when run with
        // different MPI ranks and threads.
        cell.time_sequence = arb::poisson_schedule(1e-3*params_.cell.spike_freq_hz, rng);
        return std::move(cell);
    }

    arb::cell_kind get_cell_kind(arb::cell_gid_type gid) const override {
        return arb::cell_kind::benchmark;
    }

    arb::cell_size_type num_targets(arb::cell_gid_type gid) const override {
        // Only one target, to which all incoming connections connect.
        // This could be parameterized, in which case the connections
        // generated in connections_on should end on random cell-local targets.
        return 1;
    }

    arb::cell_size_type num_sources(arb::cell_gid_type gid) const override {
        return 1;
    }

    std::vector<arb::cell_connection> connections_on(arb::cell_gid_type gid) const override {
        const auto n = params_.network.fan_in;
        std::vector<arb::cell_connection> cons;
        cons.reserve(n);
        using rng_type = std::mt19937_64;
        rng_type rng(gid);

        // Generate n incoming connections on this cell with random sources, where
        // the source can't equal gid (i.e. no self-connections).
        // We want a random distribution that will uniformly draw values from the
        // union of the two ranges: [0, gid-1] AND [gid+1, num_cells-1].
        // To do this, we draw random values in the range [0, num_cells-2], then
        // add 1 to values ≥ gid.

        std::uniform_int_distribution<arb::cell_gid_type> dist(0, params_.num_cells-2);
        for (unsigned i=0; i<n; ++i) {
            // Draw random source and adjust to avoid self-connections if neccesary.
            arb::cell_gid_type src = dist(rng);
            if (src>=gid) ++src;
            // Note: target is {gid, 0}, i.e. the first (and only) target on the cell.
            arb::cell_connection con({src, 0}, {gid, 0}, 1.f, params_.network.min_delay);
            cons.push_back(con);
        }

        return cons;
    }
};

namespace profile = arb::profile;

int main(int argc, char** argv) {
    bool is_root = true;

    try {
        arb::proc_allocation resources;
        if (auto nt = arbenv::get_env_num_threads()) {
            resources.num_threads = nt;
        }
        else {
            resources.num_threads = arbenv::thread_concurrency();
        }

#ifdef ARB_MPI_ENABLED
        arbenv::with_mpi guard(argc, argv, false);
        resources.gpu_id = arbenv::find_private_gpu(MPI_COMM_WORLD);
        auto context = arb::make_context(resources, MPI_COMM_WORLD);
        is_root = arb::rank(context) == 0;
#else
        resources.gpu_id = arbenv::default_gpu();
        auto context = arb::make_context(resources);
#endif
#ifdef ARB_PROFILE_ENABLED
        profile::profiler_initialize(context);
#endif

        std::cout << sup::mask_stream(is_root);

        bench_params params = read_options(argc, argv);

        std::cout << params << "\n";

        profile::meter_manager meters;
        meters.start(context);

        // Create an instance of our recipe.
        bench_recipe recipe(params);
        meters.checkpoint("recipe-build", context);

        // Make the domain decomposition for the model
        auto decomp = arb::partition_load_balance(recipe, context);
        meters.checkpoint("domain-decomp", context);

        // Construct the model.
        arb::simulation sim(recipe, decomp, context);
        meters.checkpoint("model-build", context);

        // Run the simulation for 100 ms, with time steps of 0.01 ms.
        sim.run(params.duration, 0.01);
        meters.checkpoint("model-run", context);

        // write meters
        auto report = profile::make_meter_report(meters, context);
        std::cout << report << "\n";

        if (is_root) {
            std::ofstream fid;
            fid.exceptions(std::ios_base::badbit | std::ios_base::failbit);
            fid.open("meters.json");
            fid << std::setw(1) << sup::to_json(report) << "\n";
        }

        // output profile and diagnostic feedback
        auto summary = profile::profiler_summary();
        std::cout << summary << "\n";

        std::cout << "there were " << sim.num_spikes() << " spikes\n";
    }
    catch (std::exception& e) {
        std::cerr << "exception caught running benchmark miniapp:\n" << e.what() << std::endl;
    }
}

std::ostream& operator<<(std::ostream& o, const bench_params& p) {
    o << "benchmark parameters:\n"
      << "  name:          " << p.name << "\n"
      << "  num cells:     " << p.num_cells << "\n"
      << "  duration:      " << p.duration << " ms\n"
      << "  fan in:        " << p.network.fan_in << " connections/cell\n"
      << "  min delay:     " << p.network.min_delay << " ms\n"
      << "  spike freq:    " << p.cell.spike_freq_hz << " Hz\n"
      << "  cell overhead: " << p.cell.realtime_ratio << " ms to advance 1 ms\n";
    o << "expected:\n"
      << "  cell advance: " << p.expected_advance_time() << " s\n"
      << "  spikes:       " << p.expected_spikes() << "\n"
      << "  events:       " << p.expected_events() << "\n"
      << "  spikes:       " << p.expected_spikes_per_interval() << " per interval\n"
      << "  events:       " << p.expected_events_per_interval()/p.num_cells << " per cell per interval";
    return o;
}

bench_params read_options(int argc, char** argv) {
    using sup::param_from_json;

    bench_params params;
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
    param_from_json(params.num_cells, "num-cells", json);
    param_from_json(params.duration, "duration", json);
    param_from_json(params.network.min_delay, "min-delay", json);
    param_from_json(params.network.fan_in, "fan-in", json);
    param_from_json(params.cell.realtime_ratio, "realtime-ratio", json);
    param_from_json(params.cell.spike_freq_hz, "spike-frequency", json);

    for (auto it=json.begin(); it!=json.end(); ++it) {
        std::cout << "  Warning: unused input parameter: \"" << it.key() << "\"\n";
    }
    std::cout << "\n";

    return params;
}
