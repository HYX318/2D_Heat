# Unit Tests CMakeLists.txt

# Include test_common.hpp directory
include_directories(${CMAKE_SOURCE_DIR}/tests)

# Collect unit test sources
set(UNIT_TEST_SOURCES
    test_array2d.cpp
    test_cartesian_topology.cpp
    test_cg_solver.cpp
    test_ghost_cell_exchange.cpp
    test_jacobi_solver.cpp
    test_mesh2d.cpp
    test_mpi_context.cpp
    test_reduction_ops.cpp
    test_sor_solver.cpp
)

# Create executable for each test or one combined test executable
add_executable(heat_equation_unit_tests ${UNIT_TEST_SOURCES})

target_link_libraries(heat_equation_unit_tests
    PRIVATE
        heat_equation_utils
        heat_equation_mpi
        heat_equation_core
        heat_equation_mesh
        ${MPI_CXX_LIBRARIES}
        gtest
        gtest_main
)

target_include_directories(heat_equation_unit_tests
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_BINARY_DIR}/include
)

# Register tests
add_test(NAME UnitTests COMMAND heat_equation_unit_tests)
