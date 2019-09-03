#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <arbor/load_balance.hpp>
#include <arbor/cable_cell.hpp>
#include <arbor/morph/morphology.hpp>
#include <arbor/morph/sample_tree.hpp>
#include <arbor/swcio.hpp>
#include <arbor/simulation.hpp>
#include <arbor/simple_sampler.hpp>

#include <sup/tinyopt.hpp>

struct options {
    std::string swc_file;
    double t_end = 20;
    double dt = 0.025;
    float syn_weight = 0.01;
};

options parse_options(int argc, char** argv);
arb::morphology default_morphology();
arb::morphology read_swc(const std::string& path);

struct single_recipe: public arb::recipe {
    explicit single_recipe(arb::morphology m): morpho(std::move(m)) {
        gprop.default_parameters = arb::neuron_parameter_defaults;
    }

    arb::cell_size_type num_cells() const override { return 1; }
    arb::cell_size_type num_probes(arb::cell_gid_type) const override { return 1; }
    arb::cell_size_type num_targets(arb::cell_gid_type) const override { return 1; }

    arb::probe_info get_probe(arb::cell_member_type probe_id) const override {
        arb::segment_location mid_soma = {0, 0.5};
        arb::cell_probe_address probe = {mid_soma, arb::cell_probe_address::membrane_voltage};

        // Probe info consists of: the probe id, a tag value to distinguish this probe
        // from others for any attached sampler (unused), and the cell probe address.

        return {probe_id, 0, probe};
    }

    arb::cell_kind get_cell_kind(arb::cell_gid_type) const override {
        return arb::cell_kind::cable;
    }

    arb::util::any get_global_properties(arb::cell_kind) const override {
        return gprop;
    }

    arb::util::unique_any get_cell_description(arb::cell_gid_type) const override {
        arb::cable_cell c = make_cable_cell(morpho, false);

        // Add HH mechanism to soma, passive channels to dendrites.
        // Discretize dendrites according to the NEURON d-lambda rule.

        c.soma()->add_mechanism("hh");
        auto& segments = c.segments();
        for (std::size_t i = 1; i<segments.size(); ++i) {
            double dx = c.segment_length_constant(100., i, gprop.default_parameters)*0.3;

            arb::cable_segment* seg = c.cable(i);
            seg->add_mechanism("pas");

            unsigned n = std::ceil(seg->length()/dx);
            seg->set_compartments(n);
        }

        // Add synapse to last branch.

        arb::cell_lid_type last_segment = morpho.num_branches()-1;
        arb::segment_location end_last_segment = { last_segment, 1. };
        c.add_synapse(end_last_segment, "exp2syn");

        return c;
    }

    arb::morphology morpho;
    arb::cable_cell_global_properties gprop;
};

int main(int argc, char** argv) {
    try {
        options opt = parse_options(argc, argv);
        single_recipe R(opt.swc_file.empty()? default_morphology(): read_swc(opt.swc_file));

        auto context = arb::make_context();
        arb::simulation sim(R, arb::partition_load_balance(R, context), context);

        // Attach a sampler to the probe described in the recipe, sampling every 0.1 ms.

        arb::trace_data<double> trace;
        sim.add_sampler(arb::all_probes, arb::regular_schedule(0.1), arb::make_simple_sampler(trace));

        // Trigger the single synapse (target is gid 0, index 0) at t = 1 ms with
        // the given weight.

        arb::spike_event spike = {{0, 0}, 1., opt.syn_weight};
        sim.inject_events({spike});

        sim.run(opt.t_end, opt.dt);

        std::cout << std::fixed << std::setprecision(4);
        for (auto entry: trace) {
            std::cout << entry.t << ", " << entry.v << "\n";
        }
    }
    catch (std::exception& e) {
        std::cerr << "caught exception: " << e.what() << "\n";
        return 2;
    }
}

options parse_options(int argc, char** argv) {
    using namespace to;
    options opt;

    char** arg = argv+1;
    while (*arg) {
        if (auto dt = parse_opt<double>(arg, 'd', "dt")) {
            opt.dt = dt.value();
        }
        else if (auto t_end = parse_opt<double>(arg, 't', "t-end")) {
            opt.t_end = t_end.value();
        }
        else if (auto weight = parse_opt<float>(arg, 'w', "weight")) {
            opt.syn_weight = weight.value();
        }
        else if (auto swc = parse_opt<std::string>(arg, 'm', "morphology")) {
            opt.swc_file = swc.value();
        }
        else {
            usage(argv[0], "[-m|--morphology SWCFILE] [-d|--dt TIME] [-t|--t-end TIME] [-w|--weight WEIGHT]");
            std::exit(1);
        }
    }
    return opt;
}

// If no SWC file is given, the default morphology consists
// of a soma of radius 6.3 µm and a single unbranched dendrite
// of length 200 µm and radius decreasing linearly from 0.5 µm
// to 0.2 µm.

arb::morphology default_morphology() {
    arb::sample_tree samples;

    auto p = samples.append(arb::msample{{  0.0, 0.0, 0.0, 6.3}, 1});
    std::cout << "p is " << p << "\n";
    p = samples.append(p, arb::msample{{  6.3, 0.0, 0.0, 0.5}, 3});
    p = samples.append(p, arb::msample{{206.3, 0.0, 0.0, 0.2}, 3});

    return arb::morphology(std::move(samples));
}

arb::morphology read_swc(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("unable to open SWC file: "+path);

    return arb::morphology(arb::swc_as_sample_tree(arb::parse_swc_file(f)));
}
