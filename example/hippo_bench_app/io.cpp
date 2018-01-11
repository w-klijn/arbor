#include <exception>
#include <fstream>
#include <iostream>
#include <istream>
#include <string>

#include <tinyopt.hpp>
#include <json/json.hpp>
#include <util/optional.hpp>
#include <util/strprintf.hpp>

#include "io.hpp"

using arb::util::optional;

namespace arb {
namespace io {


std::string usage_str = R"(
[OPTION]...

-n, --count=int        (10000)  Number of individual Poisson cell to run.

And some explanation
)";

void parse_json_options(std::string &file_name, cl_options &options);
void write_json_options(std::string &file_name, cl_options &options);

// Read options from (optional) json file and command line arguments.
cl_options read_options(int argc, char** argv, bool allow_write) {

    // The set of varuiables that might be set from the commandline
    optional<uint32_t> cells;
    optional<uint32_t> synapses_per_cell;
    optional<uint32_t> compartments_per_segment;
    optional<arb::time_type> tfinal;
    optional<std::string> json_input;
    optional<std::string> json_output;
    optional<bool> verbose;

    // Parse the possible command line parameters
    try {
        auto arg = argv + 1;
        while (*arg) {
            if (auto o = to::parse_opt<uint32_t>(arg, 'n', "cells"))                         { cells = *o; }
            else if (auto o = to::parse_opt<uint32_t>(arg, 's', "synapses_per_cell"))        { synapses_per_cell = *o; }
            else if (auto o = to::parse_opt<uint32_t>(arg, 'c', "compartments_per_segment")) { compartments_per_segment = *o; }
            else if (auto o = to::parse_opt<arb::time_type>(arg, 't', "tfinal"))             { tfinal = *o; }
            else if (auto o = to::parse_opt<bool>(arg, 'v', "verbose"))                      { verbose = *o; }
            else if (auto o = to::parse_opt<std::string>(arg, 0, "json_output"))             { json_output = *o; }
            else if (auto o = to::parse_opt<std::string>(arg, 0, "json_input"))              { json_input = *o; }
            else if (auto o = to::parse_opt(arg, 'h', "help")) {
                to::usage(argv[0], usage_str); exit(0);
            }
            else {
                throw to::parse_opt_error(*arg, "unrecognized option");
            }
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

    // Overwrite the default with 1. The json file and 2. command line options
    // Grab the default options struct
    cl_options options;

    // Read parameters from specified JSON file first, to allow
    // overriding arguments on the command line.
    if (json_input)  {
        parse_json_options(json_input.value(), options);
    }

    // Now go over the other command line arguments
    if (cells) {
        options.cells = cells.value();
    }
    if (synapses_per_cell) {
        options.synapses_per_cell = synapses_per_cell.value();
    }
    if (compartments_per_segment) {
        options.compartments_per_segment = compartments_per_segment.value();
    }
    if (tfinal) {
        options.tfinal = tfinal.value();
    }

    if (json_output && allow_write)
    {
        write_json_options(json_input.value(), options);
    }
    std::cout << "debug 4 \n";
    // If verbose output requested, emit option summary.
    if (options.verbose) {
        std::cout << options << "\n";
    }
    std::cout << "debug end io \n";
    return options;
}


void parse_json_options(std::string &file_name, cl_options &options) {
    // Ugly if statement because c++ does not have reflexion
    std::ifstream fid(file_name);
    if (fid) {
        try {
            nlohmann::json fopts;
            fid >> fopts;
            for (nlohmann::json::iterator it = fopts.begin(); it != fopts.end(); ++it) {
                // To make this if else tree readable do not follow standard code formatting
                // When adding options also add these in:  write_json_options() and  operator<<
                if (it.key() == "cells") { options.cells = it.value(); }
                else if (it.key() == "synapses_per_cell") { options.synapses_per_cell = it.value(); }
                else if (it.key() == "compartments_per_segment") { options.compartments_per_segment = it.value(); }
                else if (it.key() == "syn_type") { std::string temp = it.value(); options.syn_type = temp; }
                else if (it.key() == "morphologies") { std::string temp = it.value(); options.morphologies = temp; }
                else if (it.key() == "morph_rr") { options.morph_rr = it.value(); }
                else if (it.key() == "tfinal") { options.tfinal = it.value(); }
                else if (it.key() == "dt") { options.dt = it.value(); }
                else if (it.key() == "bin_regular") { options.bin_regular = it.value(); }
                else if (it.key() == "bin_dt") { options.bin_dt = it.value(); }
                else if (it.key() == "sample_dt") { options.sample_dt = it.value(); }
                else if (it.key() == "probe_soma_only") { options.probe_soma_only = it.value(); }
                else if (it.key() == "probe_ratio") { options.probe_ratio = it.value(); }
                else if (it.key() == "trace_prefix") { std::string temp = it.value(); options.trace_prefix = temp; }
                else if (it.key() == "trace_max_gid") { unsigned temp = it.value(); options.trace_max_gid = temp; }
                else if (it.key() == "trace_format") { std::string temp = it.value(); options.trace_format = temp; }
                else if (it.key() == "spike_file_output") { options.spike_file_output = it.value(); }
                else if (it.key() == "single_file_per_rank") { options.single_file_per_rank = it.value(); }
                else if (it.key() == "over_write") { options.over_write = it.value(); }
                else if (it.key() == "output_path") { std::string temp = it.value(); options.output_path = temp; }
                else if (it.key() == "file_name") { std::string temp = it.value(); options.file_name = temp; }
                else if (it.key() == "file_extension") { std::string temp = it.value(); options.file_extension = temp; }
                else if (it.key() == "spike_file_input") { options.spike_file_input = it.value(); }
                else if (it.key() == "input_spike_path") { std::string temp = it.value(); options.input_spike_path = temp; }
                else if (it.key() == "dry_run_ranks") { options.dry_run_ranks = it.value(); }
                else if (it.key() == "profile_only_zero") { options.profile_only_zero = it.value(); }
                else if (it.key() == "report_compartments") { options.report_compartments = it.value(); }

                else {
                    std::cerr << "Warning: Encountered an unknown key in config: " << file_name << "\n"
                        << "Key: " << it.key() << "    Value: " << it.value() << "\n";
                }
            }
        }
        catch (std::exception& e) {
            throw model_description_error(
                "unable to parse parameters in " + file_name + ": " + e.what());
        }

    }
    else {
        throw model_description_error("Unable to open file" + file_name);
    }
}

void write_json_options(std::string &file_name, cl_options &options) {
    std::ofstream fid(file_name);
    if (fid) {
        try {
            nlohmann::json fopts;
            // To make this if else tree readable do not follow standard code formatting
            fopts["cells"] = options.cells;
            fopts["synapses_per_cell"] = options.synapses_per_cell;
            fopts["compartments_per_segment"] = options.compartments_per_segment;
            fopts["syn_type"] = options.syn_type;
            if (options.morphologies) {
                fopts["morphologies"] = options.morphologies.value();
            }
            fopts["morph_rr"] = options.morph_rr;
            fopts["tfinal"] = options.tfinal;
            fopts["dt"] = options.dt;
            fopts["bin_regular"] = options.bin_regular;
            fopts["bin_dt"] = options.bin_dt;
            fopts["sample_dt"] = options.sample_dt;
            fopts["probe_soma_only"] = options.probe_soma_only;
            fopts["probe_ratio"] = options.probe_ratio;
            fopts["trace_prefix"] = options.trace_prefix;
            if (options.trace_max_gid) {
                fopts["trace_max_gid"] = options.trace_max_gid.value();
            }
            fopts["trace_format"] = options.trace_format;
            fopts["spike_file_output"] = options.spike_file_output;
            fopts["single_file_per_rank"] = options.single_file_per_rank;
            fopts["over_write"] = options.over_write;
            fopts["output_path"] = options.output_path;
            fopts["file_name"] = options.file_name;
            fopts["file_extension"] = options.file_extension;
            fopts["spike_file_input"] = options.spike_file_input;
            fopts["input_spike_path"] = options.input_spike_path;
            fopts["dry_run_ranks"] = options.dry_run_ranks;
            fopts["profile_only_zero"] = options.profile_only_zero;
            fopts["report_compartments"] = options.report_compartments;

            fid << std::setw(3) << fopts << "\n";

        }
        catch (std::exception& e) {
            throw model_description_error(
                "unable to save parameters in " + file_name + ": " + e.what());
        }
    }
    else {
        throw usage_error("unable to write to model parameter file " + file_name);
    }

}

std::ostream& operator<<(std::ostream& o, const cl_options& options) {
    cl_options defaults;

    o << "simulation options: \n";
    o << "  cells                   : " << (options.cells == defaults.cells ? " " : " * ")                                       << options.cells << "\n";
    o << "  synapses_per_cell       : " << (options.synapses_per_cell == defaults.synapses_per_cell ? " " : " * ")               << options.synapses_per_cell << "\n";
    o << "  compartments_per_segment: " << (options.compartments_per_segment == defaults.compartments_per_segment ? " " : " * ") << options.compartments_per_segment << "\n";
    o << "  syn_type                : " << (options.syn_type == defaults.syn_type ? " " : " * ")                                 << options.syn_type << "\n";
    if (options.morphologies) {
        o << "  morphologies            : " << (options.morphologies.value() == defaults.morphologies.value() ? " " : " * ")     << options.morphologies.value() << "\n";
    }
    o << "  morph_rr                : " << (options.morph_rr == defaults.morph_rr ? " " : " * ")                                 << options.morph_rr << "\n";
    o << "  tfinal                  : " << (options.tfinal == defaults.tfinal ? " " : " * ")                                     << options.tfinal << "\n";
    o << "  dt                      : " << (options.dt == defaults.dt ? " " : " * ")                                             << options.dt << "\n";
    o << "  bin_regular             : " << (options.bin_regular == defaults.bin_regular ? " " : " * ")                           << options.bin_regular << "\n";
    o << "  bin_dt                  : " << (options.bin_dt == defaults.bin_dt ? " " : " * ")                                     << options.bin_dt << "\n";
    o << "  sample_dt               : " << (options.sample_dt == defaults.sample_dt ? " " : " * ")                               << options.sample_dt << "\n";
    o << "  probe_soma_only         : " << (options.probe_soma_only == defaults.probe_soma_only ? " " : " * ")                   << options.probe_soma_only << "\n";
    o << "  probe_ratio             : " << (options.probe_ratio == defaults.probe_ratio ? " " : " * ")                           << options.probe_ratio << "\n";
    o << "  trace_prefix            : " << (options.trace_prefix == defaults.trace_prefix ? " " : " * ")                         << options.trace_prefix << "\n";
    if (options.trace_max_gid) {
        o << "  trace_max_gid           : " << (options.trace_max_gid.value() == defaults.trace_max_gid.value() ? " " : " * ")   << options.trace_max_gid.value() << "\n";
    }
    o << "  trace_format            : " << (options.trace_format == defaults.trace_format ? " " : " * ")                         << options.trace_format << "\n";
    o << "  spike_file_output       : " << (options.spike_file_output == defaults.spike_file_output ? " " : " * ")               << options.spike_file_output << "\n";
    o << "  single_file_per_rank    : " << (options.single_file_per_rank == defaults.single_file_per_rank ? " " : " * ")         << options.single_file_per_rank << "\n";
    o << "  over_write              : " << (options.over_write == defaults.over_write ? " " : " * ")                             << options.over_write << "\n";
    o << "  output_path             : " << (options.output_path == defaults.output_path ? " " : " * ")                           << options.output_path << "\n";
    o << "  file_name               : " << (options.file_name == defaults.file_name ? " " : " * ")                               << options.file_name << "\n";
    o << "  file_extension          : " << (options.file_extension == defaults.file_extension ? " " : " * ")                     << options.file_extension << "\n";
    o << "  spike_file_input        : " << (options.spike_file_input == defaults.spike_file_input ? " " : " * ")                 << options.spike_file_input << "\n";
    o << "  input_spike_path        : " << (options.input_spike_path == defaults.input_spike_path ? " " : " * ")                 << options.input_spike_path << "\n";
    o << "  dry_run_ranks           : " << (options.dry_run_ranks == defaults.dry_run_ranks ? " " : " * ")                       << options.dry_run_ranks << "\n";
    o << "  profile_only_zero       : " << (options.profile_only_zero == defaults.profile_only_zero ? " " : " * ")               << options.profile_only_zero << "\n";
    o << "  report_compartments     : " << (options.report_compartments == defaults.report_compartments ? " " : " * ")           << options.report_compartments << "\n";

    o << " \n\n Options marked with * are different from default. \n";
    return o;
}


/// Parse spike times from a stream
/// A single spike per line, trailing whitespace is ignore
/// Throws a usage error when parsing fails
///
/// Returns a vector of time_type

std::vector<time_type> parse_spike_times_from_stream(std::ifstream & fid) {
    std::vector<time_type> times;
    std::string line;
    while (std::getline(fid, line)) {
        std::stringstream s(line);

        time_type t;
        s >> t >> std::ws;

        if (!s || s.peek() != EOF) {
            throw std::runtime_error( util::strprintf(
                    "Unable to parse spike file on line %d: \"%s\"\n",
                    times.size(), line));
        }

        times.push_back(t);
    }

    return times;
}

/// Parse spike times from a file supplied in path
/// A single spike per line, trailing white space is ignored
/// Throws a usage error when opening file or parsing fails
///
/// Returns a vector of time_type

std::vector<time_type> get_parsed_spike_times_from_path(arb::util::path path) {
    std::ifstream fid(path);
    if (!fid) {
        throw std::runtime_error(util::strprintf(
            "Unable to parse spike file: \"%s\"\n", path.c_str()));
    }

    return parse_spike_times_from_stream(fid);
}

} // namespace io
} // namespace arb
