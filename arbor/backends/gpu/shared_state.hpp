#pragma once

#include <iosfwd>
#include <unordered_map>
#include <utility>
#include <vector>

#include <arbor/fvm_types.hpp>

#include "backends/gpu/gpu_store_types.hpp"

namespace arb {
namespace gpu {

/*
 * Ion state fields correspond to NMODL ion variables, where X
 * is replaced with the name of the ion. E.g. for calcium 'ca':
 *
 *     Field   NMODL variable   Meaning
 *     -------------------------------------------------------
 *     iX_     ica              calcium ion current density
 *     eX_     eca              calcium ion channel reversal potential
 *     Xi_     cai              internal calcium concentration
 *     Xo_     cao              external calcium concentration
 */

struct ion_state {
    iarray node_index_; // Instance to CV map.
    array iX_;          // (A/m²) current density
    array eX_;          // (mV) reversal potential
    array Xi_;          // (mM) internal concentration
    array Xo_;          // (mM) external concentration

    array init_Xi_;     // (mM) area-weighted initial internal concentration
    array init_Xo_;     // (mM) area-weighted initial external concentration
    array init_eX_;     // (mM) initial reversal potential

    array charge;       // charge of ionic species (global, length 1)

    ion_state() = default;

    ion_state(
        int charge,
        const std::vector<fvm_index_type>& cv,
        const std::vector<fvm_value_type>& init_Xi,
        const std::vector<fvm_value_type>& init_Xo,
        const std::vector<fvm_value_type>& init_eX,
        unsigned align
    );

    // Set ion concentrations to weighted proportion of default concentrations.
    void init_concentration();

    // Set ionic current density to zero.
    void zero_current();

    // Zero currents, reset concentrations, and reset reversal potential from
    // initial values.
    void reset();
};

struct shared_state {
    fvm_size_type n_intdom = 0; // Number of distinct integration domains.
    fvm_size_type n_cv = 0;  // Total number of CVs.
    fvm_size_type n_gj = 0;  // Total number of GJs.

    iarray cv_to_intdom;     // Maps CV index to intdom index.
    gjarray gap_junctions;   // Stores gap_junction info.
    array time;              // Maps intdom index to integration start time [ms].
    array time_to;           // Maps intdom index to integration stop time [ms].
    array dt_intdom;         // Maps intdom index to (stop time) - (start time) [ms].
    array dt_cv;             // Maps CV index to dt [ms].
    array voltage;           // Maps CV index to membrane voltage [mV].
    array current_density;   // Maps CV index to current density [A/m²].
    array conductivity;      // Maps CV index to membrane conductivity [kS/m²].

    array init_voltage;      // Maps CV index to initial membrane voltage [mV].
    array temperature_degC;  // Maps CV to local temperature (read only) [°C].

    std::unordered_map<std::string, ion_state> ion_data;

    deliverable_event_stream deliverable_events;

    shared_state() = default;

    shared_state(
        fvm_size_type n_intdom,
        const std::vector<fvm_index_type>& cv_to_intdom_vec,
        const std::vector<fvm_gap_junction>& gj_vec,
        const std::vector<fvm_value_type>& init_membrane_potential,
        const std::vector<fvm_value_type>& temperature_K,
        unsigned align
    );

    void add_ion(
        const std::string& ion_name,
        int charge,
        const std::vector<fvm_index_type>& cv,
        const std::vector<fvm_value_type>& init_iconc,
        const std::vector<fvm_value_type>& init_econc,
        const std::vector<fvm_value_type>& init_erev);

    void zero_currents();

    void ions_init_concentration();

    // Set time_to to earliest of time+dt_step and tmax.
    void update_time_to(fvm_value_type dt_step, fvm_value_type tmax);

    // Set the per-intdom and per-compartment dt from time_to - time.
    void set_dt();

    // Update gap_junction state
    void add_gj_current();

    // Return minimum and maximum time value [ms] across cells.
    std::pair<fvm_value_type, fvm_value_type> time_bounds() const;

    // Return minimum and maximum voltage value [mV] across cells.
    // (Used for solution bounds checking.)
    std::pair<fvm_value_type, fvm_value_type> voltage_bounds() const;

    // Take samples according to marked events in a sample_event_stream.
    void take_samples(
        const sample_event_stream::state& s,
        array& sample_time,
        array& sample_value);

    void reset();
};

// For debugging only
std::ostream& operator<<(std::ostream& o, shared_state& s);

} // namespace gpu
} // namespace arb
