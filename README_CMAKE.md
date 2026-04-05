# CMake Build System for 2D Heat Equation

This document describes the CMake build system for the 2D Heat Equation project.

## Requirements

- CMake 3.16 or higher
- C++17 compatible compiler (GCC 7+, Clang 5+, Apple Clang 10+)
- MPI (OpenMPI or MPICH) - required unless ENABLE_MPI=OFF
- Google Test (optional, will be downloaded automatically if not found)

## Building the Project

### Basic Build

```bash
# Create build directory
mkdir build && cd build

# Configure and build
cmake ..
make -j$(nproc)
```

### Build Options

```bash
# Disable MPI support
cmake .. -DENABLE_MPI=OFF

# Disable tests
cmake .. -DBUILD_TESTS=OFF

# Disable examples
cmake .. -DBUILD_EXAMPLES=OFF

# Combined options
cmake .. -DENABLE_MPI=ON -DBUILD_TESTS=ON -DBUILD_EXAMPLES=ON
```

### Release Build

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

### Debug Build

```bash
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j$(nproc)
```

## Build Targets

### Libraries

- `heat_equation_utils` - Utility library (logger, array2d, timer)
- `heat_equation_mpi` - MPI context library (if ENABLE_MPI=ON)

### Executables

- `timer_example` - Timer usage example (if BUILD_EXAMPLES=ON)
- `logger_test` - Logger test program (if BUILD_EXAMPLES=ON)

### Tests

- `test_array2d` - Array2D unit tests (if BUILD_TESTS=ON)
- `test_mpi_context` - MPI context unit tests (if BUILD_TESTS=ON and ENABLE_MPI=ON)

## Running Tests

```bash
# Run all tests
cd build
ctest

# Run all tests with verbose output
ctest --verbose

# Run specific test
ctest -R Array2D
ctest -R MPI_Context

# Run tests with MPI launcher (for MPI tests)
ctest --output-on-failure
```

### Running Individual Tests

```bash
# Run Array2D tests
./bin/test_array2d

# Run MPI context tests (requires MPI launcher)
mpirun -n 2 ./bin/test_mpi_context
```

## Installation

```bash
# Install to system (default: /usr/local)
sudo make install

# Install to custom prefix
cmake .. -DCMAKE_INSTALL_PREFIX=/opt/heat_equation
make install
```

## Project Structure

```
2D_Heat/
├── CMakeLists.txt                 # Root CMake configuration
├── cmake/
│   └── config.hpp.in              # Configuration file template
├── src/
│   ├── CMakeLists.txt             # Source libraries configuration
│   ├── utils/
│   │   ├── CMakeLists.txt
│   │   ├── array2d.hpp           # Header-only 2D array class
│   │   ├── timer.hpp              # Header-only timer class
│   │   ├── logger.hpp            # Logger header
│   │   └── logger.cpp            # Logger implementation
│   └── mpi/
│       ├── CMakeLists.txt
│       ├── mpi_context.hpp       # MPI context header
│       └── mpi_context.cpp       # MPI context implementation
├── tests/
│   ├── CMakeLists.txt             # Tests configuration
│   └── unit/
│       ├── CMakeLists.txt
│       ├── test_array2d.cpp      # Array2D unit tests
│       └── test_mpi_context.cpp   # MPI context unit tests
└── build/                         # Build directory (generated)
    ├── bin/                       # Executables
    ├── lib/                       # Libraries
    └── include/                   # Generated headers (config.hpp)
```

## Platform-Specific Notes

### macOS

Install required packages using Homebrew:

```bash
brew install cmake
brew install open-mpi
```

Configure with:

```bash
cmake .. -DCMAKE_CXX_COMPILER=/usr/local/bin/mpicxx
```

### Linux

Install required packages:

```bash
# Ubuntu/Debian
sudo apt-get install cmake libopenmpi-dev

# Fedora/RHEL
sudo dnf install cmake openmpi-devel
```

### Windows

Use one of the following:
- Visual Studio 2019 or later
- MinGW-w64 with MSYS2
- WSL (Windows Subsystem for Linux)

## Troubleshooting

### MPI Not Found

```bash
# Specify MPI compiler explicitly
cmake .. -DCMAKE_CXX_COMPILER=mpicxx

# Or find MPI manually
cmake .. -DMPI_CXX_COMPILER=/path/to/mpicxx
```

### Google Test Not Found

The build system will automatically download Google Test using FetchContent if not found.

### Compilation Errors

Ensure you have a C++17 compatible compiler:

```bash
# Check compiler version
g++ --version
clang++ --version
```

## Using the Libraries in Your Project

### CMake Example

```cmake
# Find the installed package
find_package(HeatEquation REQUIRED)

# Link to libraries
add_executable(my_app main.cpp)
target_link_libraries(my_app
    PRIVATE
        heat_equation_utils
        heat_equation_mpi  # Only if using MPI
)
```

### Code Example

```cpp
#include "utils/logger.hpp"
#include "utils/array2d.hpp"
#include "mpi/mpi_context.hpp"

int main(int argc, char** argv) {
    // Initialize MPI
    MPIContext mpi(argc, argv);

    // Create logger
    Logger logger;
    logger.enable_mpi_rank(true, mpi.rank());

    if (mpi.is_root()) {
        logger.info("Starting simulation");
    }

    // Create 2D array
    utils::Array2D grid(100, 100, 0.0);

    // ... simulation code ...

    return 0;
}
```

## Developer Notes

### Adding New Tests

1. Create test file in `tests/unit/`
2. Add target to `tests/unit/CMakeLists.txt`
3. Link to appropriate libraries and gtest
4. Add test with `add_test()`

### Adding New Source Files

1. Add source to appropriate `src/*/CMakeLists.txt`
2. Update library target if needed
3. Ensure headers are in correct include directories

### Code Style

The project uses strict compiler warnings:
- `-Wall`: Enable all warnings
- `-Wextra`: Enable extra warnings
- `-Wpedantic`: Warn on non-standard code

All code should compile without warnings.

## License

See LICENSE file for details.
