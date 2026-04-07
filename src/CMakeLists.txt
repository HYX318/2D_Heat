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
# Only link top-level libraries - CMake will handle transitive dependencies
# mesh depends on mpi and utils (PUBLIC), core depends on utils (PUBLIC)
target_link_libraries(heat_equation
    heat_equation_mesh
    heat_equation_core
    heat_equation_io
    ${MPI_CXX_LIBRARIES})
