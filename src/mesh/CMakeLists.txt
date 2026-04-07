# Mesh library CMakeLists.txt

# Collect source files
file(GLOB_RECURSE MESH_SOURCES *.cpp)

# Create mesh library
add_library(heat_equation_mesh ${MESH_SOURCES})
target_include_directories(heat_equation_mesh PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(heat_equation_mesh 
    PUBLIC 
        heat_equation_mpi 
        heat_equation_utils
)

# Installation
install(TARGETS heat_equation_mesh
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
