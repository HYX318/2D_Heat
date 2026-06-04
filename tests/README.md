# Testing with Google Test

This directory contains unit tests for the 2D Heat Equation Solver project, built using the Google Test framework.

## Directory Structure

```
tests/
├── README.md                    # This file
├── test_common.hpp              # Common test utilities and macros
├── CMakeLists.txt               # Test build configuration
└── unit/                        # Unit tests
    ├── CMakeLists.txt           # Unit test configuration
    ├── test_array2d.cpp
    ├── test_mesh2d.cpp
    ├── test_jacobi_solver.cpp
    ├── test_sor_solver.cpp
    ├── test_cg_solver.cpp
    └── test_mpi_context.cpp
```

## Prerequisites

- **CMake**: 3.14 or later
- **C++ Compiler**: C++17 compatible (GCC, Clang, or MSVC)
- **MPI**: OpenMPI or MPICH
- **Google Test**: Automatically downloaded via CMake FetchContent

## Building Tests

### Using CMake (Recommended)

```bash
# From project root
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTS=ON
cmake --build build
```

## Running Tests

### Using CTest

```bash
ctest --test-dir build --output-on-failure
```

CTest runs the combined unit test binary without initializing MPI. MPI-specific tests are skipped in that mode.

### MPI Test Execution

```bash
export TMPDIR=/tmp  # Required on macOS
mpirun -n 4 ./build/bin/heat_equation_unit_tests
```

### Direct Test Binary

```bash
./build/bin/heat_equation_unit_tests
```

## Test Coverage

### Array2D Tests (`test_array2d.cpp`)

14 comprehensive tests covering:

1. **Construction** - Default constructor and dimension validation
2. **FillConstruction** - Fill constructor with initial values
3. **CopyConstruction** - Deep copy functionality
4. **MoveConstruction** - Move semantics
5. **AccessorBoundsCheck** - Boundary checking for element access
6. **Fill** - Fill method with various values
7. **CopyFrom** - Copy data from another array
8. **Norms** - Max, min, L2 norm, L-infinity norm
9. **ArithmeticOperations** - +=, -=, *= operators
10. **CopyAssignment** - Copy assignment operator
11. **MoveAssignment** - Move assignment operator
12. **DataPointer** - Raw data pointer access
13. **RowMajorLayout** - Memory layout verification
14. **EdgeCases** - 1x1, 1xN, Nx1 arrays, large arrays

### MPIContext Tests (`test_mpi_context.cpp`)

8 tests covering MPI functionality:

1. **InitializationTest** - MPI initialization
2. **RankSizeTest** - Rank and size retrieval
3. **IsRootTest** - Root process identification
4. **BarrierTest** - MPI barrier synchronization
5. **NoCopyTest** - Copy semantics disabled
6. **MoveTest** - Move semantics
7. **DestructorTest** - Proper cleanup
8. **MultipleInitializationTest** - Multiple contexts

**Note**: All MPI tests require running with `mpirun`/`mpiexec`.

## Test Utilities

The `test_common.hpp` file provides useful utilities:

### Helper Functions

- `test_utils::is_mpi_initialized()` - Check MPI initialization status
- `test_utils::get_mpi_rank()` - Get current MPI rank
- `test_utils::get_mpi_size()` - Get total MPI processes
- `test_utils::is_root()` - Check if current process is root
- `test_utils::mpi_print(message)` - Print only from root process
- `test_utils::skip_if_no_mpi()` - Skip test if MPI not initialized
- `test_utils::skip_if_mpi_procs_not(n)` - Skip if process count doesn't match
- `test_utils::skip_if_mpi_procs_less_than(n)` - Skip if not enough processes

### Helper Classes

- `test_utils::Timer` - Measure test execution time
- `test_utils::MPIGuard` - RAII MPI initialization/finalization

### Custom Macros

- `ASSERT_MPI_INITIALIZED()` - Assert MPI is initialized
- `ASSERT_MPI_PROCS(n)` - Assert specific process count
- `ASSERT_IS_ROOT()` - Assert running on root process
- `ON_ROOT(code)` - Execute code only on root
- `MEASURE_TIME()` - Start timing measurement
- `PRINT_ELAPSED_TIME()` - Print elapsed time

## Writing New Tests

### Creating a New Test File

1. Create a new `.cpp` file in `tests/unit/`
2. Include Google Test and necessary headers:

```cpp
#include <gtest/gtest.h>
#include "path/to/header.hpp"
```

3. Write tests using Google Test macros:

```cpp
TEST(YourTestSuite, YourTestName) {
    // Arrange
    YourClass obj;

    // Act
    obj.doSomething();

    // Assert
    EXPECT_EQ(obj.getResult(), expected_value);
}
```

4. Update `tests/unit/CMakeLists.txt`:

```cmake
add_executable(test_your_class test_your_class.cpp)
target_link_libraries(test_your_class
    PRIVATE
        your_library
        gtest
        gtest_main
)
add_test(NAME YourClass COMMAND test_your_class)
```

### Using Test Fixtures

For tests that share common setup/teardown:

```cpp
class YourTestFixture : public ::testing::Test {
protected:
    void SetUp() override {
        // Setup code
    }

    void TearDown() override {
        // Cleanup code
    }

    YourClass* obj_;
};

TEST_F(YourTestFixture, FirstTest) {
    // Use obj_ here
}
```

### Parameterized Tests

For testing multiple inputs:

```cpp
class YourParamTest : public ::testing::TestWithParam<int> {
    // Test with int parameter
};

TEST_P(YourParamTest, DoesWork) {
    int param = GetParam();
    // Test with param
}

INSTANTIATE_TEST_SUITE_P(
    ParamRange,
    YourParamTest,
    ::testing::Values(1, 2, 3, 4, 5)
);
```

## Continuous Integration

The test framework is designed to work with CI/CD systems:

```yaml
# Example GitHub Actions workflow
- name: Build and Test
  run: |
    mkdir build && cd build
    cmake .. -DBUILD_TESTS=ON
    cmake --build . --target build_tests
    ctest --output-on-failure
```

## Troubleshooting

### MPI Tests Fail

**Problem**: MPI tests fail with "MPI not initialized" error.

**Solution**: Run with `mpirun` or `mpiexec`:

```bash
mpirun -n 4 ./build/bin/test_mpi_context
```

### Shared Memory Error (macOS)

**Problem**: Tests fail with shared memory error on macOS.

**Solution**: Set TMPDIR:

```bash
export TMPDIR=/tmp
```

### Google Test Not Found

**Problem**: CMake cannot find Google Test.

**Solution**: The build system uses FetchContent to automatically download Google Test. Ensure you have internet access during the first build.

### Build Failures

**Problem**: Compilation errors during test build.

**Solution**: Ensure you have C++17 support:

```bash
# Check compiler version
g++ --version  # or clang++ --version
```

## Additional Resources

- [Google Test Documentation](https://google.github.io/googletest/)
- [Google Test Primer](https://google.github.io/googletest/primer.html)
- [Advanced Google Test Guide](https://google.github.io/googletest/advanced.html)

## License

Tests are part of the 2D Heat Equation Solver project and follow the same license.
