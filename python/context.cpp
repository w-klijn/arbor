#include <iostream>
#include <sstream>
#include <string>

#include <pybind11/pybind11.h>

#include <arbor/context.hpp>
#include <arbor/version.hpp>
#include <arbor/util/optional.hpp>

#include "context.hpp"
#include "conversion.hpp"
#include "error.hpp"
#include "strprintf.hpp"

#ifdef ARB_MPI_ENABLED
#include "mpi.hpp"
#endif

namespace pyarb {

std::ostream& operator<<(std::ostream& o, const context_shim& ctx) {
    auto& c = ctx.context;
    const char* gpu = arb::has_gpu(c)? "True": "False";
    const char* mpi = arb::has_mpi(c)? "True": "False";
    return
        o << "<arbor.context: num_threads " << arb::num_threads(c)
          << ", has_gpu " << gpu
          << ", has_mpi " << mpi
          << ", num_ranks " << arb::num_ranks(c)
          << ">";
}

// A Python shim that holds the information that describes an arb::proc_allocation.
struct proc_allocation_shim {
    arb::util::optional<int> gpu_id = {};
    int num_threads = 1;

    proc_allocation_shim(int threads, pybind11::object gpu) {
        set_num_threads(threads);
        set_gpu_id(gpu);
    }

    proc_allocation_shim(): proc_allocation_shim(1, pybind11::none()) {}

    // getter and setter (in order to assert when being set)
    void set_gpu_id(pybind11::object gpu) {
        gpu_id = py2optional<int>(gpu, "gpu_id must be None, or a non-negative integer", is_nonneg());
    };

    void set_num_threads(int threads) {
        pyarb::assert_throw([](int n) { return n>0; }(threads), "threads must be a positive integer");
        num_threads = threads;
    };

    arb::util::optional<int> get_gpu_id() const { return gpu_id; }
    int get_num_threads() const { return num_threads; }
    bool has_gpu() const { return bool(gpu_id); }

    // helper function to use arb::make_context(arb::proc_allocation)
    arb::proc_allocation allocation() const {
        return arb::proc_allocation(num_threads, gpu_id.value_or(-1));
    }
};

std::ostream& operator<<(std::ostream& o, const proc_allocation_shim& alloc) {
    return o << "<arbor.proc_allocation: threads " << alloc.num_threads << ", gpu_id " << alloc.gpu_id << ">";
}

void register_contexts(pybind11::module& m) {
    using namespace std::string_literals;
    using namespace pybind11::literals;

    // proc_allocation
    pybind11::class_<proc_allocation_shim> proc_allocation(m, "proc_allocation",
        "Enumerates the computational resources on a node to be used for simulation.");
    proc_allocation
        .def(pybind11::init<int, pybind11::object>(),
            "threads"_a=1, "gpu_id"_a=pybind11::none(),
            "Construct an allocation with arguments:\n"
            "  threads: The number of threads available locally for execution, 1 by default.\n"
            "  gpu_id:  The identifier of the GPU to use, None by default.\n")
        .def_property("threads", &proc_allocation_shim::get_num_threads, &proc_allocation_shim::set_num_threads,
            "The number of threads available locally for execution.")
        .def_property("gpu_id", &proc_allocation_shim::get_gpu_id, &proc_allocation_shim::set_gpu_id,
            "The identifier of the GPU to use.\n"
            "Corresponds to the integer parameter used to identify GPUs in CUDA API calls.")
        .def_property_readonly("has_gpu", &proc_allocation_shim::has_gpu,
            "Whether a GPU is being used (True/False).")
        .def("__str__",  util::to_string<proc_allocation_shim>)
        .def("__repr__", util::to_string<proc_allocation_shim>);

    // context
    pybind11::class_<context_shim> context(m, "context", "An opaque handle for the hardware resources used in a simulation.");
    context
        .def(pybind11::init<>(
            [](){return context_shim(arb::make_context());}),
            "Construct a local context with one thread, no GPU, no MPI by default.\n"
            )
        .def(pybind11::init(
            [](const proc_allocation_shim& alloc){
                return context_shim(arb::make_context(alloc.allocation())); }),
            "alloc"_a,
            "Construct a local context with argument:\n"
            "  alloc:   The computational resources to be used for the simulation.\n")
#ifdef ARB_MPI_ENABLED
        .def(pybind11::init(
            [](proc_allocation_shim alloc, pybind11::object mpi){
                const char* mpi_err_str = "mpi must be None, or an MPI communicator";
                auto a = alloc.allocation(); // unwrap the C++ resource_allocation description
                if (can_convert_to_mpi_comm(mpi)) {
                    return context_shim(arb::make_context(a, convert_to_mpi_comm(mpi)));
                }
                if (auto c = py2optional<mpi_comm_shim>(mpi, mpi_err_str)) {
                    return context_shim(arb::make_context(a, c->comm));
                }
                return context_shim(arb::make_context(a));
            }),
            "alloc"_a, "mpi"_a=pybind11::none(),
            "Construct a distributed context with arguments:\n"
            "  alloc:   The computational resources to be used for the simulation.\n"
            "  mpi:     The MPI communicator, None by default.\n")
        .def(pybind11::init(
            [](int threads, pybind11::object gpu, pybind11::object mpi){
                const char* gpu_err_str = "gpu_id must be None, or a non-negative integer";
                const char* mpi_err_str = "mpi must be None, or an MPI communicator";

                auto gpu_id = py2optional<int>(gpu, gpu_err_str, is_nonneg());
                arb::proc_allocation alloc(threads, gpu_id.value_or(-1));

                if (can_convert_to_mpi_comm(mpi)) {
                    return context_shim(arb::make_context(alloc, convert_to_mpi_comm(mpi)));
                }
                if (auto c = py2optional<mpi_comm_shim>(mpi, mpi_err_str)) {
                    return context_shim(arb::make_context(alloc, c->comm));
                }

                return context_shim(arb::make_context(alloc));
            }),
            "threads"_a=1, "gpu_id"_a=pybind11::none(), "mpi"_a=pybind11::none(),
            "Construct a distributed context with arguments:\n"
            "  threads: The number of threads available locally for execution, 1 by default.\n"
            "  gpu_id:  The identifier of the GPU to use, None by default.\n"
            "  mpi:     The MPI communicator, None by default.\n")
#else
        .def(pybind11::init(
            [](int threads, pybind11::object gpu){
                auto gpu_id = py2optional<int>(gpu, "gpu_id must be None, or a non-negative integer", is_nonneg());
                return context_shim(arb::make_context(arb::proc_allocation(threads, gpu_id.value_or(-1))));
            }),
             "threads"_a=1, "gpu_id"_a=pybind11::none(),
             "Construct a local context with arguments:\n"
             "  threads: The number of threads available locally for execution, 1 by default.\n"
             "  gpu_id:  The identifier of the GPU to use, None by default.\n")
#endif
        .def_property_readonly("has_mpi", [](const context_shim& ctx){return arb::has_mpi(ctx.context);},
            "Whether the context uses MPI for distributed communication.")
        .def_property_readonly("has_gpu", [](const context_shim& ctx){return arb::has_gpu(ctx.context);},
            "Whether the context has a GPU.")
        .def_property_readonly("threads", [](const context_shim& ctx){return arb::num_threads(ctx.context);},
            "The number of threads in the context's thread pool.")
        .def_property_readonly("ranks", [](const context_shim& ctx){return arb::num_ranks(ctx.context);},
            "The number of distributed domains (equivalent to the number of MPI ranks).")
        .def_property_readonly("rank", [](const context_shim& ctx){return arb::rank(ctx.context);},
            "The numeric id of the local domain (equivalent to MPI rank).")
        .def("__str__", util::to_string<context_shim>)
        .def("__repr__", util::to_string<context_shim>);
}

} // namespace pyarb
