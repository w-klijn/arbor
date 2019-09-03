#pragma once

#include <stdexcept>
#include <string>

#include <arbor/common_types.hpp>

// Arbor-specific exception hierarchy.

namespace arb {

// Arbor internal logic error (if these are thrown,
// there is a bug in the library.)

struct arbor_internal_error: std::logic_error {
    arbor_internal_error(const std::string& what_arg):
        std::logic_error(what_arg)
    {}
};


// Common base-class for arbor run-time errors.

struct arbor_exception: std::runtime_error {
    arbor_exception(const std::string& what_arg):
        std::runtime_error(what_arg)
    {}
};

// Recipe errors:

struct bad_cell_description: arbor_exception {
    bad_cell_description(cell_kind kind, cell_gid_type gid);
    cell_gid_type gid;
    cell_kind kind;
};

struct bad_global_property: arbor_exception {
    explicit bad_global_property(cell_kind kind);
    cell_kind kind;
};

struct bad_probe_id: arbor_exception {
    explicit bad_probe_id(cell_member_type id);
    cell_member_type probe_id;
};

struct gj_kind_mismatch: arbor_exception {
    gj_kind_mismatch(cell_gid_type gid_0, cell_gid_type gid_1);
    cell_gid_type gid_0, gid_1;
};

// Domain decomposition errors:

struct gj_unsupported_domain_decomposition: arbor_exception {
    gj_unsupported_domain_decomposition(cell_gid_type gid_0, cell_gid_type gid_1);
    cell_gid_type gid_0, gid_1;
};

// Simulation errors:

struct bad_event_time: arbor_exception {
    explicit bad_event_time(time_type event_time, time_type sim_time);
    time_type event_time;
    time_type sim_time;
};

// Mechanism catalogue errors:

struct no_such_mechanism: arbor_exception {
    explicit no_such_mechanism(const std::string& mech_name);
    std::string mech_name;
};

struct duplicate_mechanism: arbor_exception {
    explicit duplicate_mechanism(const std::string& mech_name);
    std::string mech_name;
};

struct fingerprint_mismatch: arbor_exception {
    explicit fingerprint_mismatch(const std::string& mech_name);
    std::string mech_name;
};

struct no_such_parameter: arbor_exception {
    no_such_parameter(const std::string& mech_name, const std::string& param_name);
    std::string mech_name;
    std::string param_name;
};

struct invalid_parameter_value: arbor_exception {
    invalid_parameter_value(const std::string& mech_name, const std::string& param_name, const std::string& value_str);
    invalid_parameter_value(const std::string& mech_name, const std::string& param_name, double value);
    std::string mech_name;
    std::string param_name;
    std::string value_str;
    double value;
};

struct invalid_ion_remap: arbor_exception {
    explicit invalid_ion_remap(const std::string& mech_name);
    invalid_ion_remap(const std::string& mech_name, const std::string& from_ion, const std::string& to_ion);
    std::string from_ion;
    std::string to_ion;
};

struct no_such_implementation: arbor_exception {
    explicit no_such_implementation(const std::string& mech_name);
    std::string mech_name;
};

// Run-time value bounds check:

struct range_check_failure: arbor_exception {
    explicit range_check_failure(const std::string& whatstr, double value);
    double value;
};

} // namespace arb
