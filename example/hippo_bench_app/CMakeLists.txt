set(HIPPO_BENCH_SOURCES
    hippo_bench.cpp
	hippo_recipes.cpp
	../miniapp/io.cpp
    ../miniapp/morphology_pool.cpp
    ../miniapp/trace.cpp
)

add_executable(hippo_bench.exe ${HIPPO_BENCH_SOURCES})

target_link_libraries(hippo_bench.exe LINK_PUBLIC ${ARB_LIBRARIES})
target_link_libraries(hippo_bench.exe LINK_PUBLIC ${EXTERNAL_LIBRARIES})

set_target_properties(
    hippo_bench.exe
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/example/hippo_bench_app"
)

#install(FILES parse_and_plot.py
#        DESTINATION "${CMAKE_BINARY_DIR}/example/hippo_bench_app")
