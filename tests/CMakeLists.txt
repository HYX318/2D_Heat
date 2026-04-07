# Tests CMakeLists.txt

# Google Test is handled by parent CMakeLists.txt
# This file is a placeholder for future test-specific configurations

# Add unit tests
enable_testing()
add_subdirectory(unit)

# Add performance benchmarks
add_subdirectory(performance)
