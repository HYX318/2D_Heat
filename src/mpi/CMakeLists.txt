# MPI library CMakeLists.txt

file(GLOB MPI_SOURCES
    cartesian_topology.cpp
    ghost_cell_exchange.cpp
    mpi_context.cpp
    profiler.cpp
    reduction_ops.cpp
)

file(GLOB MPI_HEADERS
    cartesian_topology.hpp
    ghost_cell_exchange.hpp
    mpi_context.hpp
    profiler.hpp
    reduction_ops.hpp
)

add_library(heat_equation_mpi ${MPI_SOURCES} ${MPI_HEADERS})

target_include_directories(heat_equation_mpi
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${MPI_CXX_INCLUDE_DIRS}
)

target_link_libraries(heat_equation_mpi
    PUBLIC
        ${MPI_CXX_LIBRARIES}
    PRIVATE
        heat_equation_utils
)

# Installation
install(TARGETS heat_equation_mpi
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
install(FILES ${MPI_HEADERS}
        DESTINATION include/mpi)
