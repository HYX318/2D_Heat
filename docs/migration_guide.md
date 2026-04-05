# Migration Guide

This guide helps you migrate from older versions of the code to the current implementation.

---

# From Legacy Code to Modern Architecture

## Overview

The legacy code (Heat.cpp, HeatUtils.cpp, Interfaces.cpp) has been refactored into a modular, modern C++17 architecture with:

- **RAII pattern**: Automatic resource management
- **Template-based parallelism**: Code reuse between serial and parallel
- **Clean interfaces**: Abstract base classes and clear APIs
- **Exception safety**: Proper error handling
- **Move semantics**: Efficient object transfers

---

# API Changes

## MPI Initialization

### Before (Legacy)

```cpp
int argc, rank, size;
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);

// ... code ...

MPI_Finalize();
```

### After (Modern)

```cpp
#include "mpi/mpi_context.hpp"

int main(int argc, char** argv) {
    // RAII initialization
    MPIContext mpi(argc, argv);

    // Access rank and size
    int rank = mpi.rank();
    int size = mpi.size();
    bool is_root = mpi.is_root();

    // ... code ...

    // MPI automatically finalized when mpi goes out of scope
    return 0;
}
```

## Cartesian Topology

### Before (Legacy)

```cpp
int dims[2], periods[2], coords[2];
int neighbor_rank[4];  // N, S, E, W

// Manually create topology
MPI_Dims_create(size, 2, dims);
periods[0] = periods[1] = 0;  // Non-periodic
MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);

// Get coordinates
MPI_Cart_coords(cart_comm, 2, coords);

// Find neighbors
MPI_Cart_shift(cart_comm, 0, 1, -1, &neighbor_rank[0], &neighbor_rank[1]);  // N, S
MPI_Cart_shift(cart_comm, 1, 1, -1, &neighbor_rank[2], &neighbor_rank[3]);  // E, W

// ... code ...

MPI_Comm_free(&cart_comm);
```

### After (Modern)

```cpp
#include "mpi/cartesian_topology.hpp"

// Automatic topology creation
CartesianTopology topology(MPI_COMM_WORLD);

// Access dimensions and coordinates
int dim_x = topology.dim_x();
int dim_y = topology.dim_y();
int my_x = topology.coord_x();
int my_y = topology.coord_y();

// Access neighbors
const NeighborInfo& neighbors = topology.neighbors();
int north = neighbors.north;
int south = neighbors.south;
int east = neighbors.east;
int west = neighbors.west;

// Check boundaries
bool on_north_boundary = topology.is_on_boundary(Direction::Y, Shift::Forward);

// ... code ...

// Topology automatically freed when topology goes out of scope
```

## Ghost Cell Exchange

### Before (Legacy)

```cpp
double* u = new double[local_ny + 2][local_nx + 2];

// Exchange north-south
if (neighbor_rank[0] != MPI_PROC_NULL) {
    MPI_Sendrecv(u[1], local_nx, MPI_DOUBLE, neighbor_rank[0], 0,
                u[local_ny], local_nx, MPI_DOUBLE, neighbor_rank[0], 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

// Exchange east-west
if (neighbor_rank[2] != MPI_PROC_NULL) {
    for (int j = 0; j < local_ny + 2; ++j) {
        MPI_Sendrecv(&u[j][local_nx], 1, MPI_DOUBLE, neighbor_rank[2], 0,
                    &u[j][0], 1, MPI_DOUBLE, neighbor_rank[2], 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

// ... code ...

delete[] u;
```

### After (Modern)

```cpp
#include "mpi/ghost_cell_exchange.hpp"

// Create ghost cell exchange
GhostCellExchange exchange(local_nx, local_ny, topology);

// Create mesh with ghost cells
Mesh2D mesh(global_nx, global_ny, lx, ly, topology);

// ... update interior ...

// Simple exchange
mesh.exchange_ghost_cells();

// Or use exchange directly
exchange.exchange(mesh.data(), MPI_COMM_WORLD);

// ... code ...

// No manual cleanup needed
```

## Solver Usage

### Before (Legacy)

```cpp
double* u = new double[ny][nx];
double* u_old = new double[ny][nx];
double Coeff = 1.0 / (1.0 + 4.0 * lambda);

// Jacobi iteration
for (int iter = 0; iter < max_iter; ++iter) {
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            u[j][i] = Coeff * (
                lambda * (
                    u_old[j+1][i] + u_old[j][i+1] +
                    u_old[j-1][i] + u_old[j][i-1]
                ) + b[j][i]
            );
        }
    }

    // Swap arrays
    std::swap(u, u_old);

    // Check convergence
    // ... manual residual calculation ...
}

delete[] u;
delete[] u_old;
```

### After (Modern)

```cpp
#include "core/solver/jacobi_solver.hpp"

// Create solver
JacobiSolver<false> solver;  // Serial version

// Configure solver
SolverParams params;
params.tolerance = 1e-6;
params.max_iterations = 10000;
params.lambda = lambda;

// Solve
utils::Array2D rhs(ny, nx);
utils::Array2D solution(ny, nx);
solver.solve(rhs, solution, params);

// Get statistics
SolverStats stats = solver.get_stats();
bool converged = stats.converged;
size_t iterations = stats.iterations;
double residual = stats.final_residual;
```

## Parameter Reading

### Before (Legacy)

```cpp
int Nx, Ny, Nsteps;
double Lx, Ly, dt;
std::ifstream paramFile("Param.in");

paramFile >> Nx >> Ny;
paramFile >> Lx >> Ly;
paramFile >> dt >> Nsteps;

if (!paramFile) {
    std::cerr << "Error reading parameters" << std::endl;
    exit(1);
}
```

### After (Modern)

```cpp
#include "utils/param_reader.hpp"

ParamReader params("Param.in");

// Read with type safety and defaults
int Nx = params.get<int>("Nx", 100);
int Ny = params.get<int>("Ny", 100);
double Lx = params.get<double>("Lx", 1.0);
double Ly = params.get<double>("Ly", 1.0);
double dt = params.get<double>("dt", 0.001);
int Nsteps = params.get<int>("Nsteps", 100);

// Throws descriptive exception if file cannot be read
```

---

# Configuration File Migration

## Old Format (Param.in)

```
100 100
1.0 1.0
0.001 100
```

## New Format (params.json)

```json
{
    "grid": {
        "nx": 100,
        "ny": 100,
        "lx": 1.0,
        "ly": 1.0
    },
    "time": {
        "dt": 0.001,
        "t_final": 0.1,
        "n_steps": 100
    },
    "solver": {
        "type": "conjugate_gradient",
        "tolerance": 1e-6,
        "max_iterations": 10000,
        "preconditioner": true
    },
    "output": {
        "directory": "./output",
        "format": "csv",
        "frequency": 10
    }
}
```

## Reading New Configuration

```cpp
#include <nlohmann/json.hpp>
#include <fstream>

int main() {
    // Read JSON configuration
    std::ifstream config_file("params.json");
    nlohmann::json config;
    config_file >> config;

    // Extract parameters
    auto grid_config = config["grid"];
    size_t nx = grid_config["nx"];
    size_t ny = grid_config["ny"];
    double lx = grid_config["lx"];
    double ly = grid_config["ly"];

    auto solver_config = config["solver"];
    std::string solver_type = solver_config["type"];
    double tolerance = solver_config["tolerance"];

    // Use parameters
    // ...

    return 0;
}
```

---

# Build System Migration

## From Makefile to CMake

### Old Makefile

```makefile
CXX = mpicxx
CXXFLAGS = -O3 -std=c++11

all: heat

heat: Heat.cpp HeatUtils.cpp Interfaces.cpp
    $(CXX) $(CXXFLAGS) -o heat Heat.cpp HeatUtils.cpp Interfaces.cpp

clean:
    rm -f heat
```

### New CMakeLists.txt

```cmake
cmake_minimum_required(VERSION 3.16)
project(2D_Heat)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find MPI
find_package(MPI REQUIRED)

# Add executable
add_executable(heat_solver
    Heat.cpp
    HeatUtils.cpp
    Interfaces.cpp
)

# Link MPI
target_link_libraries(heat_solver PRIVATE MPI::MPI_CXX)

# Install target
install(TARGETS heat_solver DESTINATION bin)
```

### Building with CMake

```bash
# Configure
mkdir build && cd build
cmake ..

# Build
make -j$(nproc)

# Install
sudo make install
```

---

# Common Migration Issues and Solutions

## Issue 1: Missing MPI Finalize

**Problem**: Legacy code forgets to call `MPI_Finalize()`

**Solution**: Use `MPIContext` for automatic cleanup

```cpp
// Old code (buggy)
int main() {
    MPI_Init(&argc, &argv);
    // ... code ...
    // Forgot MPI_Finalize() - causes MPI hang
    return 0;
}

// New code (safe)
int main(int argc, char** argv) {
    MPIContext mpi(argc, argv);
    // ... code ...
    return 0;  // MPI automatically finalized
}
```

## Issue 2: Memory Leaks

**Problem**: Manual memory management leads to leaks

**Solution**: Use RAII classes

```cpp
// Old code (leaky)
double* u = new double[ny][nx];
// ... code ...
if (error) {
    return -1;  // Memory leak!
}
delete[] u;

// New code (safe)
utils::Array2D u(ny, nx);
// ... code ...
if (error) {
    return -1;  // Array automatically deleted
}
```

## Issue 3: Incorrect Array Indexing

**Problem**: Confusion between array and grid indices

**Solution**: Use `Mesh2D` which handles indexing

```cpp
// Old code (confusing)
double* u = new double[ny+2][nx+2];  // Includes ghost cells
u[j+1][i+1] = ...;  // Confusing offset

// New code (clear)
Mesh2D mesh(nx, ny, lx, ly);
mesh(j, i) = ...;  // No offset needed
```

## Issue 4: No Error Handling

**Problem**: Legacy code uses exit() on errors

**Solution**: Use exceptions for proper error handling

```cpp
// Old code (abrupt termination)
if (file.fail()) {
    std::cerr << "Error" << std::endl;
    exit(1);  // Abrupt exit
}

// New code (exception handling)
try {
    ParamReader params("params.json");
    // ... code ...
} catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;  // Clean exit
}
```

## Issue 5: Hard-Coded Parameters

**Problem**: Parameters scattered throughout code

**Solution**: Use configuration file

```cpp
// Old code (hard-coded)
constexpr int Nx = 100;
constexpr int Ny = 100;
constexpr double Lx = 1.0;

// New code (configured)
auto config = read_config("params.json");
int Nx = config["grid"]["nx"];
int Ny = config["grid"]["ny"];
double Lx = config["grid"]["lx"];
```

## Issue 6: No Convergence Checking

**Problem**: Solver runs fixed number of iterations

**Solution**: Use proper convergence criteria

```cpp
// Old code (no convergence check)
for (int iter = 0; iter < 10000; ++iter) {
    jacobi_iteration();
}

// New code (with convergence check)
JacobiariSolver solver;
SolverParams params;
params.tolerance = 1e-6;  // Convergence criterion
params.max_iterations = 10000;
solver.solve(rhs, solution, params);

if (!solver.get_stats().converged) {
    throw std::runtime_error("Solver did not converge");
}
```

---

# Step-by-Step Migration Process

## Step 1: Set Up New Build System

```bash
# Create CMakeLists.txt
cat > CMakeLists.txt << 'EOF'
cmake_minimum_required(VERSION 3.16)
project(2D_Heat)
set(CMAKE_CXX_STANDARD 17)
find_package(MPI REQUIRED)
add_executable(heat_solver Heat.cpp)
target_link_libraries(heat_solver PRIVATE MPI::MPI_CXX)
EOF

# Configure and build
mkdir build && cd build
cmake ..
make -j$(nproc)
```

## Step 2: Replace MPI Initialization

```cpp
// Add at top of main()
#include "mpi/mpi_context.hpp"

int main(int argc, char** argv) {
    MPIContext mpi(argc, argv);

    // Replace global rank/size variables
    // int rank = ...;  -> int rank = mpi.rank();
    // int size = ...;  -> int size = mpi.size();

    // ... rest of code ...
}
```

## Step 3: Replace Array Allocation

```cpp
// Add includes
#include "utils/array2d.hpp"
#include "mesh/mesh2d.hpp"

// Replace manual arrays
// double* u = new double[ny][nx];
// delete[] u;

// With automatic arrays
utils::Array2D u(ny, nx);
```

## Step 4: Replace Solver

```cpp
// Add include
#include "core/solver/jacobi_solver.hpp"

// Replace manual iteration loop
// for (int iter = 0; iter < max_iter; ++iter) { ... }

// With solver
JacobiSolver<false> solver;
SolverParams params(1e-6, 10000, 10, true, lambda);
solver.solve(rhs, solution, params);
```

## Step 5: Add Error Handling

```cpp
// Wrap main in try-catch
int main(int argc, char** argv) {
    try {
        MPIContext mpi(argc, argv);
        // ... code ...
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
```

## Step 6: Test and Verify

```bash
# Build
cd build
make -j$(nproc)

# Run tests
ctest

# Run executable
mpirun -n 4 ./heat_solver

# Compare with legacy results
diff output_new.txt output_old.txt
```

---

# Migration Checklist

Use this checklist to ensure complete migration:

- [ ] Create CMakeLists.txt
- [ ] Replace MPI_Init/Finalize with MPIContext
- [ ] Replace manual topology with CartesianTopology
- [ ] Replace manual ghost cell exchange with GhostCellExchange
- [ ] Replace manual arrays with Array2D
- [ ] Replace manual solver with ISolver implementations
- [ ] Replace manual parameter reading with ParamReader
- [ ] Add exception handling
- [ ] Update build scripts
- [ ] Update documentation
- [ ] Test serial execution
- [ ] Test parallel execution
- [ ] Verify results match legacy code
- [ ] Profile and optimize if needed
- [ ] Update configuration files

---

# Getting Help

If you encounter issues during migration:

1. **Check Documentation**: Read [API Documentation](api.md)
2. **Review Examples**: Look at files in `examples/` directory
3. **Examine Tests**: Check `tests/unit/` for usage patterns
4. **Enable Debugging**: Build with `-DCMAKE_BUILD_TYPE=Debug`
5. **Check Compiler Errors**: Read compiler messages carefully
6. **Use Valgrind**: Check for memory leaks
7. **Ask for Help**: Contact maintainers or open an issue

---

# Best Practices After Migration

1. **Use RAII**: Always use RAII classes for resource management
2. **Check Convergence**: Always verify solver convergence
3. **Handle Exceptions**: Use try-catch for error handling
4. **Profile Code**: Use profiling tools before optimization
5. **Test Thoroughly**: Test both serial and parallel versions
6. **Document Changes**: Keep documentation up to date
7. **Version Control**: Use git for tracking changes
8. **Code Review**: Have code reviewed before merging

---

# Additional Resources

- [API Documentation](api.md) - Detailed API reference
- [Examples Guide](examples.md) - Complete code examples
- [Architecture Documentation](architecture.md) - System design
- [Performance Guide](performance_guide.md) - Optimization tips
