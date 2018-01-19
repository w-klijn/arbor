set(HIPPO_BENCH_SOURCES
    hippo_bench.cpp
    io.cpp
    hippo_bench_recipes.cpp
    morphology_pool.cpp
    trace.cpp
)

add_executable(hippo_bench.exe ${HIPPO_BENCH_SOURCES})

target_link_libraries(hippo_bench.exe LINK_PUBLIC ${ARB_LIBRARIES})
target_link_libraries(hippo_bench.exe LINK_PUBLIC ${EXTERNAL_LIBRARIES})

if(ARB_WITH_MPI)
    target_link_libraries(hippo_bench.exe LINK_PUBLIC ${MPI_C_LIBRARIES})
    set_property(TARGET hippo_bench.exe APPEND_STRING PROPERTY LINK_FLAGS "${MPI_C_LINK_FLAGS}")
endif()

set_target_properties(
    hippo_bench.exe
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/example"
)
