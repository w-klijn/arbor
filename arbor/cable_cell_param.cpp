#include <cfloat>
#include <cmath>

#include <arbor/cable_cell_param.hpp>

#include "util/maputil.hpp"

namespace arb {

void check_global_properties(const cable_cell_global_properties& G) {
    auto& param = G.default_parameters;

    if (!param.init_membrane_potential) {
        throw cable_cell_error("missing global default parameter value: init_membrane_potential");
    }

    if (!param.temperature_K) {
        throw cable_cell_error("missing global default parameter value: temperature");
    }

    if (!param.axial_resistivity) {
        throw cable_cell_error("missing global default parameter value: axial_resistivity");
    }

    if (!param.membrane_capacitance) {
        throw cable_cell_error("missing global default parameter value: membrane_capacitance");
    }

    for (const auto& ion: util::keys(G.ion_species)) {
        if (!param.ion_data.count(ion)) {
            throw cable_cell_error("missing ion defaults for ion "+ion);
        }
    }

    for (const auto& kv: param.ion_data) {
        auto& ion = kv.first;
        const cable_cell_ion_data& data = kv.second;
        if (std::isnan(data.init_int_concentration)) {
            throw cable_cell_error("missing init_int_concentration for ion "+ion);
        }
        if (std::isnan(data.init_ext_concentration)) {
            throw cable_cell_error("missing init_ext_concentration for ion "+ion);
        }
        if (std::isnan(data.init_reversal_potential) && !param.reversal_potential_method.count(ion)) {
            throw cable_cell_error("missing init_reversal_potential or reversal_potential_method for ion "+ion);
        }
    }
}

cable_cell_local_parameter_set neuron_parameter_defaults = {
    // ion defaults: internal concentration [mM], external concentration [mM], reversal potential [mV]
    {
        {"na", {10.0,  140.0,  115 - 65.}},
        {"k",  {54.4,    2.5,  -12 - 65.}},
        {"ca", {5e-5,    2.0,  12.5*std::log(2.0/5e-5)}}
    },
    // initial membrane potential [mV]
    -65.0,
    // temperatue [K]
    6.3 + 273.15,
    // axial resistivity [Ω·cm]
    35.4,
    // membrane capacitance [F/m²]
    0.01
};

} // namespace arb
