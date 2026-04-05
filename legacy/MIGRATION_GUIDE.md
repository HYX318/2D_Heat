# Migration Guide: Legacy Code to New Architecture

## Overview

This guide helps you migrate from the legacy API to the new 2D Heat Equation solver architecture.

## Why Migrate?

The new architecture provides:
- **Modern C++17 features**: RAII, smart pointers, type safety
- **Better performance**: Optimized memory layout and communication patterns
- **Easier to use**: Cleaner API with less boilerplate
- **More maintainable**: Clear separation of concerns, better code organization
- **Testable**: Unit tests and integration tests

## Quick Reference

| Legacy Function | New Equivalent | Header |
|----------------|----------------|---------|
| `ReadParam()` | Use configuration file or parameters | See examples |
| `Contiguous2D()` | `utils::Array2D` or `Mesh2D` | `src/utils/array2d.hpp` |
| `Init()` | `Mesh2D::apply_bc()` | `src/mesh/mesh2d.hpp` |
| `ExactSol()` | Custom function or `Mesh2D::fill()` | `src/mesh/mesh2d.hpp` |
| `Analytic()` | Lambda function | Inline |
| `Error()` | `Array2D::operator-=` | `src/utils/array2d.hpp` |
| `TwoNorm()` | `Array2D::l2_norm()` | `src/utils/array2d.hpp` |
| `InftyNorm()` | `Array2D::linfty_norm()` | `src/utils/array2d.hpp` |
| `Copy()` | `Array2D::copy_from()` | `src/utils/array2d.hpp` |
| `Jacobi()` | `JacobiSolver<false>` | `src/core/solver/jacobi_solver.hpp` |
| `MPIJacobi()` | `JacobiSolver<true>` | `src/core/solver/jacobi_solver.hpp` |
| `Interfaces()` | `GhostCellExchange::exchange()` | `src/mpi/ghost_cell_exchange.hpp` |
| `Export()` | Custom I/O or examples | See examples |

## Detailed Migration Examples

### 1. Array Allocation

**Legacy:**
```cpp
int Nx = 100, Ny = 100;
double** array = new double*[Ny+2];
Contiguous2D(array, Nx+2, Ny+2);
// ... use array ...
delete[] array[0];
delete[] array;
```

**New:**
```cpp
#include "src/utils/array2d.hpp"

int Nx = 100, Ny = 100;
utils::Array2D array(Ny+2, Nx+2);  // Note: rows, cols order
// ... use array ...
// Automatic cleanup (RAII)
```

**Accessing Elements:**
```cpp
// Legacy: array[j][i]
// New:   array(j, i)
```

### 2. Mesh Creation

**Legacy:**
```cpp
int Nx = 100, Ny = 100;
double* x = new double[Nx+2];
double* y = new double[Ny+2];
// Compute x, y arrays...
double** Sol = new double*[Ny+2];
Contiguous2D(Sol, Nx+2, Ny+2);
Init(Sol, x, y, Nx, Ny);
```

**New:**
```cpp
#include "src/mesh/mesh2d.hpp"

int Nx = 100, Ny = 100;
double lx = 1.0, ly = 1.0;

// Serial version
Mesh2D mesh(Nx, Ny, lx, ly);

// Parallel version
CartesianTopology topology(MPI_COMM_WORLD);
Mesh2D mesh(Nx, Ny, lx, ly, topology);

// Initialize with analytic solution
auto init_func = [](double x, double y, double t) {
    return std::sin(M_PI*x) * std::sin(2*M_PI*y);
};
mesh.apply_bc(init_func, 0.0);
```

### 3. Solver Usage

**Legacy (Serial):**
```cpp
double** x = new double*[Ny+2];
double** x0 = new double*[Ny+2];
double** b = new double*[Ny+2];
// ... allocate arrays ...
double lambda = 0.25;
double Residu;
double Tol = 1e-6;
int iConv;

Jacobi(x, x0, b, Residu, Tol, iConv, Nx, Ny, lambda);
```

**New (Serial):**
```cpp
#include "src/core/solver/jacobi_solver.hpp"

utils::Array2D x(Ny+2, Nx+2);
utils::Array2D x0(Ny+2, Nx+2);
utils::Array2D b(Ny+2, Nx+2);
// ... fill b ...

double lambda = 0.25;
SolverParams params;
params.tolerance = 1e-6;
params.max_iterations = 10000;
params.lambda = lambda;
params.verbose = true;

JacobiSolver<false> solver;  // Serial solver
solver.solve(b, x, params);

auto stats = solver.get_stats();
std::cout << "Iterations: " << stats.iterations << std::endl;
std::cout << "Residual: " << stats.final_residual << std::endl;
```

**Legacy (Parallel):**
```cpp
MPI_Comm SBD_COMM;
MPI_Datatype colType;
int NeighbourRank[4];
// ... setup MPI topology ...

MPIJacobi(x, x0, b, Residu, Tol, iConv, Nx, Ny, lambda,
          NeighbourRank, SBD_COMM, colType, myRank);
```

**New (Parallel):**
```cpp
#include "src/mpi/cartesian_topology.hpp"
#include "src/mpi/ghost_cell_exchange.hpp"
#include "src/core/solver/jacobi_solver.hpp"

// Setup topology
CartesianTopology topology(MPI_COMM_WORLD);

// Create mesh with ghost cells
Mesh2D mesh(Nx, Ny, lx, ly, topology);

// Create ghost cell exchange
GhostCellExchange ghost_exchange(mesh.nx(), mesh.ny(), topology);

// Setup solver
SolverParams params;
params.tolerance = 1e-6;
params.max_iterations = 10000;
params.lambda = lambda;
params.verbose = (topology.rank() == 0);

JacobiSolver<true> solver(&ghost_exchange, MPI_COMM_WORLD);
solver.solve(rhs, solution, params);

auto stats = solver.get_stats();
```

### 4. Error Computation

**Legacy:**
```cpp
double** Err = new double*[Ny+2];
Contiguous2D(Err, Nx+2, Ny+2);
Error(Sol, Ref, Err, Nx, Ny);

double h = 1.0/(Nx+1);
double l2_err = TwoNorm(Err, Nx, Ny, h);
double linfty_err = InftyNorm(Err, Nx, Ny);
```

**New:**
```cpp
utils::Array2D Err = Sol;
Err -= Ref;  // Element-wise subtraction

// Note: TwoNorm in legacy includes h factor
double h = 1.0/(Nx+1);
double l2_err = Err.l2_norm() * h;
double linfty_err = Err.linfty_norm();
```

### 5. MPI Communication

**Legacy:**
```cpp
MPI_Comm SBD_COMM;
MPI_Datatype colType;
int NeighbourRank[4];
// ... setup ...

Interfaces(Sol, NeighbourRank, SBD_COMM, colType, myRank, Nx, Ny);
```

**New:**
```cpp
CartesianTopology topology(MPI_COMM_WORLD);
Mesh2D mesh(Nx, Ny, lx, ly, topology);

// Automatic ghost cell exchange
mesh.exchange_ghost_cells();
```

## Step-by-Step Migration

### Step 1: Replace Array Allocation

1. Find all `Contiguous2D` calls
2. Replace with `utils::Array2D`
3. Update index access from `[j][i]` to `(j, i)`
4. Remove manual memory deallocation

### Step 2: Replace Mesh Setup

1. Find all `x` and `y` coordinate array allocations
2. Replace with `Mesh2D` class
3. Replace `Init()` with `apply_bc()`
4. Use mesh coordinate methods instead of separate arrays

### Step 3: Replace Solvers

1. Find `Jacobi` calls → Replace with `JacobiSolver<false>`
2. Find `MPIJacobi` calls → Replace with `JacobiSolver<true>`
3. Create `SolverParams` structure
4. Update error checking to use `get_stats()`

### Step 4: Replace MPI Setup

1. Find manual `MPI_Cart_create` calls
2. Replace with `CartesianTopology` class
3. Replace `Interfaces()` calls with `GhostCellExchange`
4. Use mesh's built-in ghost cell support

### Step 5: Replace Error Computation

1. Replace `Error()` with array arithmetic
2. Replace `TwoNorm()` with `l2_norm()`
3. Replace `InftyNorm()` with `linfty_norm()`
4. Adjust for any scaling factors

### Step 6: Remove Legacy Headers

1. Remove `#include "HeatUtils.h"`
2. Remove `#include` "Interfaces.h"`
3. Remove `#include "Consts.h"`
4. Add new architecture headers

### Step 7: Update CMakeLists.txt

```cmake
# Old
target_sources(your_target
    PRIVATE
    HeatUtils.cpp
    Interfaces.cpp
)

# New
target_sources(your_target
    PRIVATE
    src/utils/array2d.cpp
    src/mesh/mesh2d.cpp
    src/mpi/cartesian_topology.cpp
    src/mpi/ghost_cell_exchange.cpp
    src/core/solver/jacobi_solver.cpp
)
```

## Common Issues and Solutions

### Issue 1: Array Index Order

**Problem:** Legacy uses `array[j][i]` but new uses `array(i, j)`

**Solution:** Remember that legacy's first index is row (j), second is column (i).
The new API uses (row, col) which is consistent.

### Issue 2: Ghost Cell Offsets

**Problem:** Legacy manually handles ghost cell indices (0 and N+1)

**Solution:** New `Mesh2D` class handles this automatically:
- Interior points: `mesh(i, j)` where i=0..nx-1, j=0..ny-1
- Ghost cells: `mesh.at(j, i)` for direct access

### Issue 3: MPI Datatype Creation

**Problem:** Legacy manually creates `colType` for MPI

**Solution:** `GhostCellExchange` handles this automatically:
```cpp
GhostCellExchange ghost_exchange(nx, ny, topology);
MPI_Datatype colType = ghost_exchange.col_type();  // If needed
```

### Issue 4: Neighbor Information

**Problem:** Legacy uses `NeighbourRank[4]` array

**Solution:** Use `CartesianTopology::neighbors()`:
```cpp
const auto& neighbors = topology.neighbors();
int south_rank = neighbors.south;
int north_rank = neighbors.north;
// ...
```

## Complete Example: Legacy to New

### Legacy Version
```cpp
// Heat.cpp - simplified
#include "HeatUtils.h"
#include "Interfaces.h"
#include "Consts.h"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    // Read parameters
    int Nx, Ny, Nt;
    double StabP;
    ReadParam(Nx, Ny, Nt, StabP);

    // Create Cartesian topology
    MPI_Comm SBD_COMM;
    int dims[2] = {3, 3};
    int periods[2] = {0, 0};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &SBD_COMM);

    // Get neighbors
    int NeighbourRank[4];
    MPI_Cart_shift(SBD_COMM, 0, 1, &NeighbourRank[West], &NeighbourRank[East]);
    MPI_Cart_shift(SBD_COMM, 1, 1, &NeighbourRank[South], &NeighbourRank[North]);

    // Allocate arrays
    double** Sol = new double*[Ny+2];
    Contiguous2D(Sol, Nx+2, Ny+2);
    double** Sol0 = new double*[Ny+2];
    Contiguous2D(Sol0, Nx+2, Ny+2);
    double** b = new double*[Ny+2];
    Contiguous2D(b, Nx+2, Ny+2);

    // Initialize
    double* x = new double[Nx+2];
    double* y = new double[Ny+2];
    // ... compute x, y ...
    Init(Sol, x, y, Nx, Ny);

    // Solve
    double Residu, lambda = 0.25, Tol = 1e-6;
    int iConv;
    Copy(Sol, b, Nx, Ny);
    Copy(Sol, Sol0, Nx, Ny);

    MPI_Datatype colType;
    MPI_Type_vector(Ny+2, 1, Nx+2, MPI_DOUBLE, &colType);
    MPI_Type_commit(&colType);

    MPIJacobi(Sol, Sol0, b, Residu, Tol, iConv, Nx, Ny, lambda,
              NeighbourRank, SBD_COMM, colType, 0);

    // Cleanup
    delete[] Sol[0]; delete[] Sol;
    delete[] Sol0[0]; delete[] Sol0;
    delete[] b[0]; delete[] b;
    delete[] x; delete[] y;
    MPI_Type_free(&colType);
    MPI_Finalize();

    return 0;
}
```

### New Version
```cpp
// new_main.cpp - equivalent
#include "src/mesh/mesh2d.hpp"
#include "src/mpi/cartesian_topology.hpp"
#include "src/mpi/ghost_cell_exchange.hpp"
#include "src/core/solver/jacobi_solver.hpp"
#include "src/utils/logger.hpp"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    // Read parameters (use new configuration system)
    int Nx = 100, Ny = 100, Nt = 100;
    double StabP = 0.25;

    // Create topology
    CartesianTopology topology(MPI_COMM_WORLD);

    // Create mesh with ghost cells
    double lx = 1.0, ly = 1.0;
    Mesh2D mesh(Nx, Ny, lx, ly, topology);

    // Initialize with analytic solution
    auto init_func = [](double x, double y, double t) {
        return std::sin(M_PI*x) * std::sin(2*M_PI*y);
    };
    mesh.apply_bc(init_func, 0.0);

    // Create ghost cell exchange
    GhostCellExchange ghost_exchange(Nx, Ny, topology);

    // Setup solver
    utils::Array2D rhs = mesh.data();
    utils::Array2D solution(Ny+2, Nx+2);

    SolverParams params;
    params.tolerance = 1e-6;
    params.max_iterations = 10000;
    params.lambda = StabP;
    params.verbose = (topology.rank() == 0);

    // Solve
    JacobiSolver<true> solver(&ghost_exchange, MPI_COMM_WORLD);
    solver.solve(rhs, solution, params);

    auto stats = solver.get_stats();
    if (topology.rank() == 0) {
        LOG_INFO("Converged in " << stats.iterations << " iterations");
        LOG_INFO("Final residual: " << stats.final_residual);
    }

    MPI_Finalize();
    return 0;
}
```

## Performance Considerations

The new architecture provides several performance benefits:

1. **Memory Management**: RAII prevents leaks and improves cache locality
2. **Communication**: `GhostCellExchange` uses optimized non-blocking communication
3. **Solver**: Template-based solvers enable compiler optimizations
4. **Vectorization**: Contiguous memory layout enables SIMD optimizations

## Getting Help

- See examples in `examples/` directory
- Check unit tests in `tests/unit/` directory
- Review API documentation in header files
- Check implementation summaries in docs/

## Deprecation Timeline

- **Version 1.0.0**: Legacy compatibility layer available (current)
- **Version 1.1.0**: Legacy functions will produce warnings at compile time
- **Version 2.0.0**: Legacy compatibility layer will be removed

We recommend migrating to the new API as soon as possible.
