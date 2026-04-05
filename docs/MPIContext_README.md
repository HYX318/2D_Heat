# MPIContext Class Documentation

## Overview

`MPIContext` is a RAII (Resource Acquisition) wrapper for the Message Passing Interface (MPI) library. It automatically handles MPI initialization and finalization, preventing common errors and providing a clean, exception-safe interface for MPI operations.

## Features

- **RAII Pattern**: Automatic initialization on construction and cleanup on destruction
- **Thread Safety**: Supports MPI threading with `MPI_THREAD_FUNNELED`
- **Move-Only Semantics**: Prevents accidental copies, allows moves for ownership transfer
- **Exception Safety**: Handles MPI errors gracefully with descriptive messages
- **Double Initialization Handling**: Correctly manages cases where MPI is already initialized

## File Structure

```
src/mpi/
├── mpi_context.hpp    # Header file with class declaration
└── mpi_context.cpp    # Implementation file

tests/unit/
└── test_mpi_context.cpp  # Unit tests
```

## Usage Example

### Basic Usage

```cpp
#include "mpi/mpi_context.hpp"

int main(int argc, char** argv) {
    // Initialize MPI (RAII - automatic cleanup)
    MPIContext mpi(argc, argv);

    // Get process information
    int rank = mpi.rank();
    int size = mpi.size();

    if (mpi.is_root()) {
        std::cout << "Running with " << size << " processes" << std::endl;
    }

    // Synchronize all processes
    mpi.barrier();

    // MPI is automatically finalized when mpi goes out of scope
    return 0;
}
```

### Error Handling

```cpp
try {
    MPIContext mpi(argc, argv);
    // ... MPI operations
} catch (const std::runtime_error& e) {
    std::cerr << "MPI Error: " << e.what() << std::endl;
    return 1;
}
```

### Move Semantics

```cpp
// Move constructor
MPIContext context1(argc, argv);
MPIContext context2(std::move(context1));

// Move assignment
MPIContext context3(argc, argv);
context3 = std::move(context2);
```

## Class Interface

### Constructor

```cpp
MPIContext(int& argc, char** argv)
```

Initializes MPI if not already initialized using `MPI_Init_thread` with `MPI_THREAD_FUNNELED` support. Retrieves process rank and size.

**Throws**: `std::runtime_error` if MPI initialization fails.

### Destructor

```cpp
~MPIContext()
```

Automatically finalizes MPI if this context was responsible for initialization and MPI hasn't been finalized yet.

### Accessor Methods

```cpp
int rank() const          // Returns current process rank (0 to size-1)
int size() const          // Returns total number of processes
bool is_root() const     // Returns true if rank == 0
```

### Control Methods

```cpp
void barrier() const
// Calls MPI_Barrier to synchronize all processes
// Throws: std::runtime_error if barrier fails or context is finalized

void abort(int error_code, const std::string& message) const
// Calls MPI_Abort to terminate all MPI processes
// Root process prints the message before aborting
```

## Compilation

### Building the class

```bash
mpic++ -Wall -g -c src/mpi/mpi_context.cpp -o mpi_context.o
```

### Building tests

```bash
# Use the provided build script
./scripts/build_mpi_context_test.sh
```

Or manually:

```bash
mkdir -p build
mpic++ -Wall -g -c src/mpi/mpi_context.cpp -o build/mpi_context.o
mpic++ -Wall -g -std=c++11 build/mpi_context.o tests/unit/test_mpi_context.cpp \
  -o build/test_mpi_context -lgtest -lgtest_main -lpthread
```

## Running Tests

### Using the provided script

```bash
# Run with 4 processes (default)
./scripts/run_mpi_context_test.sh

# Run with specified number of processes
./scripts/run_mpi_context_test.sh 8
```

### Using mpirun directly

```bash
# Set TMPDIR for macOS to avoid shared memory errors
export TMPDIR=/tmp

# Run tests
mpirun -np 4 ./build/test_mpi_context
```

## Test Coverage

The unit tests cover the following aspects:

1. **InitializationTest**: Verifies MPI initialization and context creation
2. **RankSizeTest**: Validates rank and size retrieval
3. **IsRootTest**: Tests root process identification
4. **BarrierTest**: Tests barrier synchronization
5. **NoCopyTest**: Verifies copy semantics are disabled
6. **MoveTest**: Tests move semantics and ownership transfer
7. **DestructorTest**: Tests automatic cleanup
8. **MultipleInitializationTest**: Tests handling of existing MPI initialization

## Requirements

- MPI library (OpenMPI, MPICH, or compatible)
- C++11 or later
- Google Test framework for unit tests

## macOS-Specific Notes

On macOS, you may encounter shared memory errors with MPI. Set the TMPDIR environment variable:

```bash
export TMPDIR=/tmp
```

Add this to your `.zshrc` or `.bashrc` for persistence.

## Design Considerations

### Why RAII?

- Prevents resource leaks (MPI_Finalize not called)
- Makes ownership explicit
- Exception-safe construction and destruction
- Fits naturally with C++ scope-based resource management

### Why Move-Only?

- Only one context should own MPI initialization
- Prevents accidental duplication that could lead to double finalization
- Allows transfer of ownership when needed

### Thread Support

Uses `MPI_THREAD_FUNNELED` which:
- Allows only the main thread to make MPI calls
- Provides the best performance for typical MPI applications
- Is sufficient for most parallel algorithms

## Best Practices

1. Create `MPIContext` early in `main()`
2. Store by reference when passing to functions
3. Use `is_root()` to conditionally execute I/O operations
4. Call `barrier()` when timing critical sections
5. Let RAII handle cleanup - don't call `MPI_Finalize()` manually

## Common Pitfalls

- **Don't create multiple contexts**: While the class handles it, design for single ownership
- **Don't use after move**: Moved-from objects are in invalid state
- **Don't copy**: Copy semantics are explicitly disabled
- **Don't call MPI functions directly**: Use the context methods when available

## License

This code is part of the 2D Heat Equation Solver project.
