# Performance Tests / Benchmarks CMakeLists.txt

# Collect performance benchmark sources
set(PERFORMANCE_BENCHMARK_SOURCES
    benchmark_communication.cpp
    benchmark_mpi_communication.cpp
    benchmark_overlap.cpp
    benchmark_scaling.cpp
)

# Create executable for benchmarks
add_executable(heat_equation_benchmarks ${PERFORMANCE_BENCHMARK_SOURCES})

target_link_libraries(heat_equation_benchmarks
    PRIVATE
        heat_equation_utils
        heat_equation_mpi
        heat_equation_core
        heat_equation_mesh
        ${MPI_CXX_LIBRARIES}
)

target_include_directories(heat_equation_benchmarks
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_BINARY_DIR}/include
)

# Register benchmark as a test (can be run separately)
add_test(NAME PerformanceBenchmarks COMMAND heat_equation_benchmarks)
