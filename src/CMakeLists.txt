# Main source CMakeLists.txt

# Collect all source files
set(SOURCES
    main.cpp
)

# Utilities library (must come first as other modules depend on it)
add_subdirectory(utils)

# MPI library
add_subdirectory(mpi)

# Core library (solvers and schemes)
add_subdirectory(core)

# I/O library
add_subdirectory(io)

# Mesh library
add_subdirectory(mesh)

# Build main executable
add_executable(heat_equation ${SOURCES})
target_link_libraries(heat_equation
    heat_equation_options
    heat_equation_utils
    heat_equation_mesh
    heat_equation_mpi
    heat_equation_core
    heat_equation_io
    ${MPI_CXX_LIBRARIES})
