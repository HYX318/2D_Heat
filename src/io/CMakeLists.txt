# I/O library CMakeLists.txt

# Collect source files
file(GLOB_RECURSE IO_SOURCES *.cpp)

# Create I/O library
add_library(heat_equation_io ${IO_SOURCES})
target_include_directories(heat_equation_io 
    PUBLIC 
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${MPI_CXX_INCLUDE_DIRS}
)
target_link_libraries(heat_equation_io 
    PRIVATE 
        heat_equation_utils
        heat_equation_mesh
        ${MPI_CXX_LIBRARIES}
)
