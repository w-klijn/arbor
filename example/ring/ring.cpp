/*
 * A miniapp that demonstrates how to make a ring model
 *
 */

#include <fstream>
#include <iomanip>
#include <iostream>

#include <nlohmann/json.hpp>

#include <arbor/assert_macro.hpp>
#include <arbor/common_types.hpp>
#include <arbor/cable_cell.hpp>
#include <arbor/context.hpp>
#include <arbor/load_balance.hpp>
#include <arbor/cable_cell.hpp>
#include <arbor/profile/meter_manager.hpp>
#include <arbor/profile/profiler.hpp>
#include <arbor/simple_sampler.hpp>
#include <arbor/simulation.hpp>
#include <arbor/recipe.hpp>
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

// Parameters used to generate the random cell morphologies.
struct cell_parameters {
    cell_parameters() = default;

    //  Maximum number of levels in the cell (not including the soma)
    unsigned max_depth = 5;

    // The following parameters are described as ranges.
    // The first value is at the soma, and the last value is used on the last level.
    // Values at levels in between are found by linear interpolation.
    std::array<double,2> branch_probs = {1.0, 0.5}; //  Probability of a branch occuring.
    std::array<unsigned,2> compartments = {20, 2};  //  Compartment count on a branch.
    std::array<double,2> lengths = {200, 20};       //  Length of branch in μm.

    // The number of synapses per cell.
    unsigned synapses = 1;
};

struct ring_params {
    ring_params() = default;

    std::string name = "default";
    unsigned num_cells = 10;
    double min_delay = 10;
    double duration = 100;
    cell_parameters cell;
};

ring_params read_options(int argc, char** argv);
using arb::cell_gid_type;
using arb::cell_lid_type;
using arb::cell_size_type;
using arb::cell_member_type;
using arb::cell_kind;
using arb::time_type;
using arb::cell_probe_address;

// Writes voltage trace as a json file.
void write_trace_json(const arb::trace_data<double>& trace);

// Generate a cell.
arb::cable_cell branch_cell(arb::cell_gid_type gid, const cell_parameters& params);

class ring_recipe: public arb::recipe {
public:
    ring_recipe(unsigned num_cells, cell_parameters params, unsigned min_delay):
        num_cells_(num_cells),
        cell_params_(params),
        min_delay_(min_delay)
    {
        gprop_.default_parameters = arb::neuron_parameter_defaults;
    }

    cell_size_type num_cells() const override {
        return num_cells_;
    }

    arb::util::unique_any get_cell_description(cell_gid_type gid) const override {
        return branch_cell(gid, cell_params_);
    }

    cell_kind get_cell_kind(cell_gid_type gid) const override {
        return cell_kind::cable;
    }

    // Each cell has one spike detector (at the soma).
    cell_size_type num_sources(cell_gid_type gid) const override {
        return 1;
    }

    // The cell has one target synapse, which will be connected to cell gid-1.
    cell_size_type num_targets(cell_gid_type gid) const override {
        return 1;
    }

    // Each cell has one incoming connection, from cell with gid-1.
    std::vector<arb::cell_connection> connections_on(cell_gid_type gid) const override {
        cell_gid_type src = gid? gid-1: num_cells_-1;
        return {arb::cell_connection({src, 0}, {gid, 0}, event_weight_, min_delay_)};
    }

    // Return one event generator on gid 0. This generates a single event that will
    // kick start the spiking.
    std::vector<arb::event_generator> event_generators(cell_gid_type gid) const override {
        std::vector<arb::event_generator> gens;
        if (!gid) {
            gens.push_back(arb::explicit_generator(arb::pse_vector{{{0, 0}, 1.0, event_weight_}}));
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

    arb::util::any get_global_properties(arb::cell_kind) const override {
        return gprop_;
    }


private:
    cell_size_type num_cells_;
    cell_parameters cell_params_;
    double min_delay_;
    float event_weight_ = 0.05;
    arb::cable_cell_global_properties gprop_;
};

struct cell_stats {
    using size_type = unsigned;
    size_type ncells = 0;
    size_type nsegs = 0;
    size_type ncomp = 0;

    cell_stats(arb::recipe& r) {
#ifdef ARB_MPI_ENABLED
        int nranks, rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nranks);
        ncells = r.num_cells();
        size_type cells_per_rank = ncells/nranks;
        size_type b = rank*cells_per_rank;
        size_type e = (rank==nranks-1)? ncells: (rank+1)*cells_per_rank;
        size_type nsegs_tmp = 0;
        size_type ncomp_tmp = 0;
        for (size_type i=b; i<e; ++i) {
            auto c = arb::util::any_cast<arb::cable_cell>(r.get_cell_description(i));
            nsegs_tmp += c.num_segments();
            ncomp_tmp += c.num_compartments();
        }
        MPI_Allreduce(&nsegs_tmp, &nsegs, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&ncomp_tmp, &ncomp, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
#else
        ncells = r.num_cells();
        for (size_type i=0; i<ncells; ++i) {
            auto c = arb::util::any_cast<arb::cable_cell>(r.get_cell_description(i));
            nsegs += c.num_segments();
            ncomp += c.num_compartments();
        }
#endif
    }

    friend std::ostream& operator<<(std::ostream& o, const cell_stats& s) {
        return o << "cell stats: "
                 << s.ncells << " cells; "
                 << s.nsegs << " segments; "
                 << s.ncomp << " compartments.";
    }
};

int main(int argc, char** argv) {
    try {
        bool root = true;

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
        root = arb::rank(context) == 0;
#else
        resources.gpu_id = arbenv::default_gpu();
        auto context = arb::make_context(resources);
#endif

#ifdef ARB_PROFILE_ENABLED
        arb::profile::profiler_initialize(context);
#endif

        std::cout << sup::mask_stream(root);

        // Print a banner with information about hardware configuration
        std::cout << "gpu:      " << (has_gpu(context)? "yes": "no") << "\n";
        std::cout << "threads:  " << num_threads(context) << "\n";
        std::cout << "mpi:      " << (has_mpi(context)? "yes": "no") << "\n";
        std::cout << "ranks:    " << num_ranks(context) << "\n" << std::endl;

        auto params = read_options(argc, argv);

        arb::profile::meter_manager meters;
        meters.start(context);

        // Create an instance of our recipe.
        ring_recipe recipe(params.num_cells, params.cell, params.min_delay);
        cell_stats stats(recipe);
        std::cout << stats << "\n";

        auto decomp = arb::partition_load_balance(recipe, context);

        // Construct the model.
        arb::simulation sim(recipe, decomp, context);

        // Set up the probe that will measure voltage in the cell.

        // The id of the only probe on the cell: the cell_member type points to (cell 0, probe 0)
        auto probe_id = cell_member_type{0, 0};
        // The schedule for sampling is 10 samples every 1 ms.
        auto sched = arb::regular_schedule(0.1);
        // This is where the voltage samples will be stored as (time, value) pairs
        arb::trace_data<double> voltage;
        // Now attach the sampler at probe_id, with sampling schedule sched, writing to voltage
        sim.add_sampler(arb::one_probe(probe_id), sched, arb::make_simple_sampler(voltage));

        // Set up recording of spikes to a vector on the root process.
        std::vector<arb::spike> recorded_spikes;
        if (root) {
            sim.set_global_spike_callback(
                [&recorded_spikes](const std::vector<arb::spike>& spikes) {
                    recorded_spikes.insert(recorded_spikes.end(), spikes.begin(), spikes.end());
                });
        }

        meters.checkpoint("model-init", context);

        std::cout << "running simulation" << std::endl;
        // Run the simulation for 100 ms, with time steps of 0.025 ms.
        sim.run(params.duration, 0.025);

        meters.checkpoint("model-run", context);

        auto ns = sim.num_spikes();

        // Write spikes to file
        if (root) {
            std::cout << "\n" << ns << " spikes generated at rate of "
                      << params.duration/ns << " ms between spikes\n";
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

        // Write the samples to a json file.
        if (root) write_trace_json(voltage);

        auto report = arb::profile::make_meter_report(meters, context);
        std::cout << report;
    }
    catch (std::exception& e) {
        std::cerr << "exception caught in ring miniapp: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

void write_trace_json(const arb::trace_data<double>& trace) {
    std::string path = "./voltages.json";

    nlohmann::json json;
    json["name"] = "ring demo";
    json["units"] = "mV";
    json["cell"] = "0.0";
    json["probe"] = "0";

    auto& jt = json["data"]["time"];
    auto& jy = json["data"]["voltage"];

    for (const auto& sample: trace) {
        jt.push_back(sample.t);
        jy.push_back(sample.v);
    }

    std::ofstream file(path);
    file << std::setw(1) << json << "\n";
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

    // Add additional synapses that will not be connected to anything.
    for (unsigned i=1u; i<params.synapses; ++i) {
        cell.add_synapse({1, 0.5}, "expsyn");
    }

    return cell;
}

ring_params read_options(int argc, char** argv) {
    using sup::param_from_json;

    ring_params params;
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
    param_from_json(params.min_delay, "min-delay", json);
    param_from_json(params.cell.max_depth, "depth", json);
    param_from_json(params.cell.branch_probs, "branch-probs", json);
    param_from_json(params.cell.compartments, "compartments", json);
    param_from_json(params.cell.lengths, "lengths", json);
    param_from_json(params.cell.synapses, "synapses", json);

    if (!json.empty()) {
        for (auto it=json.begin(); it!=json.end(); ++it) {
            std::cout << "  Warning: unused input parameter: \"" << it.key() << "\"\n";
        }
        std::cout << "\n";
    }

    return params;
}

