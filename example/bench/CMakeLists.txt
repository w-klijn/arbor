add_executable(bench EXCLUDE_FROM_ALL bench.cpp)
add_dependencies(examples bench)

target_link_libraries(bench PRIVATE arbor arborenv arbor-sup ext-tclap ext-json)
