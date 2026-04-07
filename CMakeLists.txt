cmake_minimum_required(VERSION 3.16)
project(HeatEquation VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
endif()

# Build options
option(BUILD_TESTS "Build unit tests" OFF)
option(BUILD_EXAMPLES "Build example programs" ON)
option(ENABLE_MPI "Enable MPI support" ON)

# MPI support
if(ENABLE_MPI)
    find_package(MPI REQUIRED)
    message(STATUS "MPI FOUND: ${MPI_FOUND}")
    message(STATUS "MPI COMPILER: ${MPI_CXX_COMPILER}")
    message(STATUS "MPI VERSION: ${MPI_VERSION}")
else()
    message(STATUS "MPI support disabled")
endif()

# Google Test (only if tests are enabled)
if(BUILD_TESTS)
    # Try to find GTest first - much faster than downloading
    find_package(GTest QUIET)

    if(NOT GTest_FOUND)
        message(STATUS "GTest not found, will use FetchContent to download it")
        message(WARNING "GTest download can take several minutes - requires internet connection")
        include(FetchContent)
        FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG v1.14.0
        )
        # For Windows: Prevent overriding of parent project's compiler/linker settings
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(googletest)
    else()
        message(STATUS "GTest found: ${GTEST_VERSION}")
        message(STATUS "GTest library: ${GTEST_LIBRARIES}")
    endif()
endif()

# Generate config.hpp
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/include/config.hpp
)

# Add subdirectories
add_subdirectory(src)
# add_subdirectory(legacy)  # Legacy compatibility layer (disabled for now)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Installation rules
install(TARGETS heat_equation_utils heat_equation_mesh heat_equation_mpi
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)

install(DIRECTORY src/utils/ src/mesh/ src/mpi/
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp")

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/config.hpp
        DESTINATION include)
