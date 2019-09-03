#include <vector>

#include <arbor/mechanism.hpp>
#include <arbor/version.hpp>

#include "backends/multicore/fvm.hpp"

#ifdef ARB_GPU_ENABLED
#include "backends/gpu/fvm.hpp"
#endif

#include "common.hpp"
#include "mech_private_field_access.hpp"
#include "fvm_lowered_cell.hpp"
#include "fvm_lowered_cell_impl.hpp"
#include "sampler_map.hpp"
#include "simple_recipes.hpp"
#include "unit_test_catalogue.hpp"

using namespace arb;

using backend = arb::multicore::backend;
using fvm_cell = arb::fvm_lowered_cell_impl<backend>;

using shared_state = backend::shared_state;
ACCESS_BIND(std::unique_ptr<shared_state> fvm_cell::*, private_state_ptr, &fvm_cell::state_)

template <typename backend>
void run_test(std::string mech_name,
        std::vector<std::string> state_variables,
        std::unordered_map<std::string, fvm_value_type> assigned_variables,
        std::vector<fvm_value_type> t0_values,
        std::vector<fvm_value_type> t1_values) {

    auto cat = make_unit_test_catalogue();

    fvm_size_type ncell = 1;
    fvm_size_type ncv = 1;
    std::vector<fvm_index_type> cv_to_intdom(ncv, 0);

    std::vector<fvm_gap_junction> gj = {};
    auto instance = cat.instance<backend>(mech_name);
    auto& test = instance.mech;

    std::vector<fvm_value_type> temp(ncv, 300.);
    std::vector<fvm_value_type> vinit(ncv, -65);

    auto shared_state = std::make_unique<typename backend::shared_state>(
            ncell, cv_to_intdom, gj, vinit, temp, test->data_alignment());

    mechanism_layout layout;
    mechanism_overrides overrides;

    layout.weight.assign(ncv, 1.);
    for (fvm_size_type i = 0; i<ncv; ++i) {
        layout.cv.push_back(i);
    }

    test->instantiate(0, *shared_state, overrides, layout);

    for (auto a: assigned_variables) {
        test->set_parameter(a.first, std::vector<fvm_value_type>(ncv,a.second));
    }

    shared_state->reset();

    test->initialize();

    if (!t0_values.empty()) {
        for (unsigned i = 0; i < state_variables.size(); i++) {
            for (unsigned j = 0; j < ncv; j++) {
                EXPECT_NEAR(t0_values[i], mechanism_field(test.get(), state_variables[i]).at(j), 1e-6);
            }
        }
    }

    shared_state->update_time_to(0.5, 0.5);
    shared_state->set_dt();

    test->nrn_state();

    if (!t1_values.empty()) {
        for (unsigned i = 0; i < state_variables.size(); i++) {
            for (unsigned j = 0; j < ncv; j++) {
                EXPECT_NEAR(t1_values[i], mechanism_field(test.get(), state_variables[i]).at(j), 1e-6);
            }
        }
    }
}

TEST(mech_kinetic, kintetic_scaled) {
    std::vector<std::string> state_variables = {"s", "h", "d"};
    std::vector<fvm_value_type> t0_values = {0.5, 0.2, 0.3};
    std::vector<fvm_value_type> t1_0_values = {0.373297, 0.591621, 0.0350817};
    std::vector<fvm_value_type> t1_1_values = {0.329897, 0.537371, 0.132732};

    run_test<multicore::backend>("test0_kin_compartment", state_variables, {}, t0_values, t1_0_values);
    run_test<multicore::backend>("test1_kin_compartment", state_variables, {}, t0_values, t1_1_values);

}

TEST(mech_kinetic, kintetic_1_conserve) {
    std::vector<std::string> state_variables = {"s", "h", "d"};
    std::vector<fvm_value_type> t0_values = {0.5, 0.2, 0.3};
    std::vector<fvm_value_type> t1_values = {0.380338, 0.446414, 0.173247};

    run_test<multicore::backend>("test0_kin_diff", state_variables, {}, t0_values, t1_values);
    run_test<multicore::backend>("test0_kin_conserve", state_variables, {}, t0_values, t1_values);
}

TEST(mech_kinetic, kintetic_2_conserve) {
    std::vector<std::string> state_variables = {"a", "b", "x", "y"};
    std::vector<fvm_value_type> t0_values = {0.2, 0.8, 0.6, 0.4};
    std::vector<fvm_value_type> t1_values = {0.217391304, 0.782608696, 0.33333333, 0.66666666};

    run_test<multicore::backend>("test1_kin_diff", state_variables, {}, t0_values, t1_values);
    run_test<multicore::backend>("test1_kin_conserve", state_variables, {}, t0_values, t1_values);
}

TEST(mech_linear, linear) {
    std::vector<std::string> state_variables = {"h", "s", "d"};
    std::vector<fvm_value_type> values = {0.5, 0.2, 0.3};
    std::unordered_map<std::string, fvm_value_type> assigned_variables = {{"a0", 2.5}, {"a1",0.5}, {"a2",3}, {"a3",2.3}};

    run_test<multicore::backend>("test_linear_state", state_variables, assigned_variables, {}, values);
    run_test<multicore::backend>("test_linear_init", state_variables, assigned_variables, values, {});
    run_test<multicore::backend>("test_linear_init_shuffle", state_variables, assigned_variables, values, {});
}

#ifdef ARB_GPU_ENABLED
TEST(mech_kinetic_gpu, kintetic_scaled) {
     std::vector<std::string> state_variables = {"s", "h", "d"};
    std::vector<fvm_value_type> t0_values = {0.5, 0.2, 0.3};
    std::vector<fvm_value_type> t1_0_values = {0.373297, 0.591621, 0.0350817};
    std::vector<fvm_value_type> t1_1_values = {0.329897, 0.537371, 0.132732};

    run_test<gpu::backend>("test0_kin_compartment", state_variables, {}, t0_values, t1_0_values);
    run_test<gpu::backend>("test1_kin_compartment", state_variables, {}, t0_values, t1_1_values);
}

TEST(mech_kinetic_gpu, kintetic_1_conserve) {
    std::vector<std::string> state_variables = {"s", "h", "d"};
    std::vector<fvm_value_type> t0_values = {0.5, 0.2, 0.3};
    std::vector<fvm_value_type> t1_values = {0.380338, 0.446414, 0.173247};

    run_test<gpu::backend>("test0_kin_diff", state_variables, {}, t0_values, t1_values);
    run_test<gpu::backend>("test0_kin_conserve", state_variables, {}, t0_values, t1_values);
}

TEST(mech_kinetic_gpu, kintetic_2_conserve) {
    std::vector<std::string> state_variables = {"a", "b", "x", "y"};
    std::vector<fvm_value_type> t0_values = {0.2, 0.8, 0.6, 0.4};
    std::vector<fvm_value_type> t1_values = {0.217391304, 0.782608696, 0.33333333, 0.66666666};

    run_test<gpu::backend>("test1_kin_diff", state_variables, {}, t0_values, t1_values);
    run_test<gpu::backend>("test1_kin_conserve", state_variables, {}, t0_values, t1_values);
}

TEST(mech_linear_gpu, linear) {
    std::vector<std::string> state_variables = {"h", "s", "d"};
    std::vector<fvm_value_type> values = {0.5, 0.2, 0.3};
    std::unordered_map<std::string, fvm_value_type> assigned_variables = {{"a0", 2.5},{"a1",0.5},{"a2",3},{"a3",2.3}};

    run_test<gpu::backend>("test_linear_state", state_variables, assigned_variables, {}, values);
    run_test<gpu::backend>("test_linear_init", state_variables, assigned_variables, values, {});
    run_test<gpu::backend>("test_linear_init_shuffle", state_variables, assigned_variables, values, {});
}

#endif
