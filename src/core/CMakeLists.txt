# Core library CMakeLists.txt

# Collect all source files
file(GLOB_RECURSE SOLVER_SOURCES solver/*.cpp)
file(GLOB_RECURSE SCHEME_SOURCES scheme/*.cpp)

# Create core library
add_library(heat_equation_core ${SOLVER_SOURCES} ${SCHEME_SOURCES})
target_include_directories(heat_equation_core 
    PUBLIC 
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${MPI_CXX_INCLUDE_DIRS}
)
target_link_libraries(heat_equation_core 
    PUBLIC 
        heat_equation_utils
        ${MPI_CXX_LIBRARIES}
)
