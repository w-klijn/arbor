#pragma once

#include <string>
#include <vector>

#include "backends/event.hpp"
#include "backends/multicore/matrix_state.hpp"
#include "backends/multicore/multi_event_stream.hpp"
#include "backends/multicore/multicore_common.hpp"
#include "backends/multicore/shared_state.hpp"
#include "backends/multicore/threshold_watcher.hpp"
#include "execution_context.hpp"
#include "util/padded_alloc.hpp"
#include "util/range.hpp"
#include "util/rangeutil.hpp"

namespace arb {
namespace multicore {

struct backend {
    static bool is_supported() { return true; }
    static std::string name() { return "cpu"; }

    using value_type = fvm_value_type;
    using index_type = fvm_index_type;
    using size_type  = fvm_size_type;

    using array  = arb::multicore::array;
    using iarray = arb::multicore::iarray;

    static util::range<const value_type*> host_view(const array& v) {
        return util::range_pointer_view(v);
    }

    static util::range<const index_type*> host_view(const iarray& v) {
        return util::range_pointer_view(v);
    }

    using matrix_state = arb::multicore::matrix_state<value_type, index_type>;
    using threshold_watcher = arb::multicore::threshold_watcher;

    using deliverable_event_stream = arb::multicore::deliverable_event_stream;
    using sample_event_stream = arb::multicore::sample_event_stream;

    using shared_state = arb::multicore::shared_state;

    static threshold_watcher voltage_watcher(
        const shared_state& state,
        const std::vector<index_type>& cv,
        const std::vector<value_type>& thresholds,
        const execution_context& context)
    {
        return threshold_watcher(
            state.cv_to_intdom.data(),
            state.time.data(),
            state.time_to.data(),
            state.voltage.data(),
            cv,
            thresholds,
            context);
    }
};

} // namespace multicore
} // namespace arb
