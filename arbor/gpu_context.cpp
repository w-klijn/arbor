#include <memory>

#include <arbor/arbexcept.hpp>

#include "gpu_context.hpp"

#ifdef ARB_HAVE_GPU
#include <cuda.h>
#include <cuda_runtime.h>
#endif

namespace arb {

enum gpu_flags {
    has_concurrent_managed_access = 1,
    has_atomic_double = 2
};

gpu_context_handle make_gpu_context(int id) {
    return std::make_shared<gpu_context>(id);
}

bool gpu_context_has_gpu(const gpu_context& ctx) {
    return ctx.has_gpu();
}

bool gpu_context::has_concurrent_managed_access() const {
    return attributes_ & gpu_flags::has_concurrent_managed_access;
}

bool gpu_context::has_atomic_double() const {
    return attributes_ & gpu_flags::has_atomic_double;
}

bool gpu_context::has_gpu() const {
    return id_ != -1;
}

#ifndef ARB_HAVE_GPU

void gpu_context::set_gpu() const {
    throw arbor_exception("Arbor must be compiled with CUDA support to set a GPU.");
}

void gpu_context::synchronize_for_managed_access() const {}

gpu_context::gpu_context(int) {
    throw arbor_exception("Arbor must be compiled with CUDA support to select a GPU.");
}

#else

gpu_context::gpu_context(int gpu_id) {
    cudaDeviceProp prop;
    auto status = cudaGetDeviceProperties(&prop, gpu_id);
    if (status==cudaErrorInvalidDevice) {
        throw arbor_exception("Invalid GPU id " + std::to_string(gpu_id));
    }

    // Set the device.
    // The device will also have to be set for every host thread that uses the
    // GPU, however performing this call here is a good check that the GPU can
    // be set and initialized.
    id_ = gpu_id;
    set_gpu();

    // Record the device attributes
    attributes_ = 0;
    if (prop.concurrentManagedAccess) {
        attributes_ |= gpu_flags::has_concurrent_managed_access;
    }
    if (prop.major*100 + prop.minor >= 600) {
        attributes_ |= gpu_flags::has_atomic_double;
    }
}

void gpu_context::synchronize_for_managed_access() const {
    if(!has_concurrent_managed_access()) {
        cudaDeviceSynchronize();
    }
}

void gpu_context::set_gpu() const {
    if (!has_gpu()) {
        throw arbor_exception(
            "Call to gpu_context::set_gpu() when there is no GPU selected.");
    }
    auto status = cudaSetDevice(id_);
    if (status != cudaSuccess) {
        throw arbor_exception(
            "Unable to select GPU id " + std::to_string(id_)
            + ": " + cudaGetErrorName(status));
    }
}

#endif

} // namespace arb
