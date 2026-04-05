# Quick Start Guide

## Prerequisites

1. **CMake** (3.16+)
   ```bash
   # macOS
   brew install cmake

   # Ubuntu/Debian
   sudo apt-getting install cmake

   # Fedora/RHEL
   sudo dnf install cmake
   ```

2. **MPI** (OpenMPI or MPICH)
   ```bash
   # macOS
   brew install open-mpi

   # Ubuntu/Debian
   sudo apt-get install libopenmpi-dev

   # Fedora/RHEL
   sudo dnf install openmpi-devel
   ```

3. **C++17 Compiler** (GCC 7+, Clang 5+, Apple)
   - Usually already installed with your system

## Building

### 1. Configure and Build

```bash
# Create build directory
mkdir build && cd build

# Configure (Release mode with optimizations)
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build
make -j$(nproc)
```

### 2. Run Tests

```bash
# Run all tests
ctest

# Run tests with verbose output
ctest --verbose

# Run specific test
ctest -R Array2D
```

### 3. Build Options

```bash
# Disable MPI (single-threaded build)
cmake .. -DENABLE_MPI=OFF

# Disable tests
cmake .. -DBUILD_TESTS=OFF

# Debug build
cmake .. -DCMAKE_BUILD_TYPE=Debug
```

## Using the Libraries

### Example Code

```cpp
#include "utils/logger.hpp"
#include "utils/array2d.hpp"
#include "utils/timer.hpp"
#include "mpi/mpi_context.hpp"

int main(int argc, char** argv) {
    // Initialize MPI
    MPIContext mpi(argc, argv);

    // Create logger with MPI rank
    Logger logger;
    logger.enable_mpi_rank(true, mpi.rank());

    if (mpi.is_root()) {
        logger.info("Starting 2D heat equation simulation");
    }

    // Create 2D grid
    utils::Array2D grid(100, 100, 0.0);

    // Timer
    Timer timer;
    timer.start();

    // Simulation code here...

    timer.stop();
    if (mpi.is_root()) {
        logger.info("Simulation completed in " +
                   std::to_string(timer.elapsed_ms()) + " ms");
    }

    return 0;
}
```

### CMake Integration

```cmake
cmake_minimum_required(VERSION 3.16)
project(MyHeatApp)

set(CMAKE_CXX_STANDARD 17)

# Find MPI
find_package(MPI REQUIRED)

# Add executable
add_executable(my_app main.cpp)

# Link to project libraries
target_link_libraries(my_app
    PRIVATE
        heat_equation_utils
        heat_equation_mpi
        MPI::MPI_CXX
)

target_include_directories(my_app
    PRIVATE
        /path/to/2D_Heat/src
)
```

## Project Structure

```
2D_Heat/
├── src/
│   ├── utils/           # Utility classes
│   │   ├── array2d.hpp   # 2D array with mathematical operations
│   │   ├── timer.hpp      # High-precision timer
│   │   ├── logger.hpp     # Threaded logger
│   │   └── logger.cpp
│   └── mpi/             # MPI utilities
│       └── mpi_context.hpp  # RAII MPI wrapper
├── tests/              # Unit tests
│   └── unit/
│       ├── test_array2d.cpp
│       └── test_mpi_context.cpp
└── build/              # Build directory (generated)
    ├── bin/             # Executables
    └── lib/             # Libraries
```

## Common Issues

### MPI Not Found

```bash
# Specify MPI compiler explicitly
cmake .. -DCMAKE_CXX_COMPILER=mpicxx
```

### Build Errors

```bash
# Clean and rebuild
rm -rf build
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j$(nproc)
```

### Tests Not Running

```bash
# Make sure tests were built
cmake .. -DBUILD_TESTS=ON
make -j$(nproc)

# Run tests
ctest --output-on-failure
```

## Next Steps

- Read `README_CMAKE.md` for detailed CMake documentation
- Check header files in `src/utils/` and `src/mpi/` for API details
- Look at test files in `tests/unit/` for usage examples
- Run examples in `build/bin/` directory

## Support

For issues or questions:
1. Check the README files
2. Review the test cases
3. Examine the header file documentation
