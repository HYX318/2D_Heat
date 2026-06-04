# Performance Tests / Benchmarks CMakeLists.txt

# Collect performance benchmark sources - each has its own main(), build separately
# benchmark_communication.cpp
add_executable(benchmark_communication benchmark_communication.cpp)
target_link_libraries(benchmark_communication
    PRIVATE
        heat_equation_utils
        heat_equation_mpi
        heat_equation_core
        heat_equation_mesh
        ${MPI_CXX_LIBRARIES}
)
target_include_directories(benchmark_communication
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_BINARY_DIR}/include
)

# benchmark_mpi_communication.cpp
add_executable(benchmark_mpi_communication benchmark_mpi_communication.cpp)
target_link_libraries(benchmark_mpi_communication
    PRIVATE
        heat_equation_utils
        heat_equation_mpi
        heat_equation_core
        heat_equation_mesh
        ${MPI_CXX_LIBRARIES}
)
target_include_directories(benchmark_mpi_communication
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_BINARY_DIR}/include
)

# benchmark_overlap.cpp
add_executable(benchmark_overlap benchmark_overlap.cpp)
target_link_libraries(benchmark_overlap
    PRIVATE
        heat_equation_utils
        heat_equation_mpi
        heat_equation_core
        heat_equation_mesh
        ${MPI_CXX_LIBRARIES}
)
target_include_directories(benchmark_overlap
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_BINARY_DIR}/include
)

# benchmark_scaling.cpp
add_executable(benchmark_scaling benchmark_scaling.cpp)
target_link_libraries(benchmark_scaling
    PRIVATE
        heat_equation_utils
        heat_equation_mpi
        heat_equation_core
        heat_equation_mesh
        ${MPI_CXX_LIBRARIES}
)
target_include_directories(benchmark_scaling
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_BINARY_DIR}/include
)
