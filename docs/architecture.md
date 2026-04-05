# Architecture Documentation

## Project Overview

The 2D Heat Equation Solver is a high-performance parallel computing application designed to solve the 2D heat diffusion equation using MPI (Message Passing Interface). The project employs modern C++17 design patterns with a focus on modularity, performance, and maintainability.

### Mathematical Model

The solver implements the heat equation:

```
∂u/∂t = α (∂²u/∂x² + ∂²u/∂y²)
```

Discretized using implicit time integration schemes:
- Implicit Euler: (I - λ∇²) u^{n+1} = u^n
- Crank-Nicolson: (I - λ/2∇²) u^{n+1} = (I + λ/2∇²) u^n

Where λ = αΔt/Δh² is the diffusion coefficient.

## Architecture Design Principles

### 1. Layered Architecture

The codebase follows a strict layered architecture:

```
┌─────────────────────────────────────┐
│     Application Layer              │  (main.cpp, examples)
├─────────────────────────────────────┤
│     Time Stepping Layer            │  (ImplicitEuler, CrankNicolson)
├─────────────────────────────────────┤
│     Solver Layer                   │  (ISolver, Jacobi, SOR, CG)
├─────────────────────────────────────┤
│     Mesh Layer                     │  (Mesh2D, DomainDecomposition)
├─────────────────────────────────────┤
│     MPI Layer                      │  (MPIContext, CartesianTopology)
├─────────────────────────────────────┤
│     Utilities Layer                │  (Array2D, Logger, Timer)
└─────────────────────────────────────┘
```

### 2. RAII Pattern

All resources follow the RAII (Resource Acquisition Is Initialization) pattern:
- Automatic initialization on construction
- Automatic cleanup on destruction
- No manual resource management required

### 3. Move Semantics

Classes implement move semantics for efficient transfers:
- Move constructors for fast object transfer
- Move assignment operators
- Copy operations deleted where appropriate

### 4. Exception Safety

The codebase is exception-safe:
- Strong guarantee on critical operations
- No resource leaks on exceptions
- Clear error messages

### 5. Template-Based Parallelism

Parallel and serial versions share the same implementation:
- Template parameter controls MPI communication
- Code reuse between serial and parallel versions
- Compile-time optimization

## Module Description

### MPI Layer

**Location**: `src/mpi/`

The MPI layer provides all parallel computing infrastructure:

#### MPIContext
- RAII wrapper for MPI initialization/finalization
- Automatic cleanup on destruction
- Thread-safe initialization support
- Provides rank, size, and barrier operations

#### CartesianTopology
- Manages 2D Cartesian topology for domain decomposition
- Automatic dimension calculation (closest to square)
- Neighbor identification for communication patterns
- Boundary detection

#### GhostCellExchange
- Efficient ghost cell exchange between processes
- Non-blocking communication support
- Automatic boundary handling
- Overlapping communication with computation

#### ReductionOps
- Type-safe MPI reduction operations
- Support for scalar and array reductions
- Special operations: scan, exscan, maxloc, minloc
- Broadcast support

#### Profiler
- Performance profiling for MPI operations
- Timing and statistics collection
- Communication pattern analysis

### Mesh Layer

**Location**: `src/mesh/`

The mesh layer handles grid management and domain decomposition:

#### Mesh2D
- 2D rectangular mesh with ghost cell support
- Dirichlet boundary conditions
- Coordinate transformations (local/global, grid/physical)
- Numerical operations (Laplacian, norms, arithmetic)
- Efficient 5-point stencil operations

Grid layout with ghost cells:
```
Row 0      : South ghost row
Row [1..ny]: Interior rows
Row ny+1   : North ghost row

Col 0      : West ghost column
Col [1..nx]: Interior columns
Col nx+1   : East ghost column
```

#### DomainDecomposition
- Divides global domain into subdomains
- Supports automatic and manual decomposition
- Coordinate transformation between global and local
- Load balance analysis

#### CoordinateSystem
- Uniform and non-uniform grid generation
- Coordinate indexing and lookup
- Coordinate transformations
- Grid quality metrics

### Solver Layer

**Location**: `src/core/solver/`

The solver layer implements various iterative methods:

#### ISolver (Interface)
- Abstract base class for all solvers
- Common interface: solve(), get_stats(), reset()
- Solver parameter and statistics structures
- Enables runtime solver selection

#### JacobiSolver
- Jacobi iterative method implementation
- Template-based serial/parallel versions
- Efficient ghost cell exchange
- Residual computation with global reduction

#### SORSolver
- Successive Over-Relaxation method
- Automatic omega optimization
- Red-black ordering for parallelism
- Gauss-Seidel special case (omega=1.0)

#### ConjugateGradientSolver
- Conjugate Gradient method
- Optional Jacobi preconditioning
- Restart mechanism for numerical stability
- Superior convergence for Poisson-like problems

### Time Stepping Layer

**Location**: `src/core/timestepping/`

The time stepping layer implements time integration schemes:

#### ImplicitEuler
- First-order accurate in time
- Unconditionally stable
- Solves (I - λ∇²) u^{n+1} = u^n
- Simple and robust

#### CrankNicolson
- Second-order accurate in time
- Unconditionally stable
- Solves (I - λ/2∇²) u^{n+1} = (I + λ/2∇²) u^n
- Better accuracy for given timestep

### Utilities Layer

**Location**: `src/utils/`

The utilities layer provides common functionality:

#### Array2D
- RAII 2D array wrapper
- Row-major contiguous memory layout
- Bounds-checked access
- Mathematical operations (norms, arithmetic)

#### Logger
- Thread-safe structured logging
- Multiple log levels (DEBUG, INFO, WARNING, ERROR, FATAL)
- Timestamp and MPI rank support
- Multiple output targets (console, file)

#### Timer
- High-precision timing
- Start/stop functionality
- Elapsed time queries

#### ParamReader
- Configuration file reader
- Type-safe parameter access
- Default value support

#### SolutionExporter
- Result file writer
- Various output formats (CSV, binary)
- Automatic file naming

## Data Flow Diagram

```
┌─────────────┐
│ Input Params│
└──────┬──────┘
       │
       ▼
┌─────────────┐     ┌─────────────┐
│   CMake     │────▶│  MPI Init   │
│  Config     │     └──────┬──────┘
└─────────────┘            │
                          ▼
┌─────────────┐     ┌─────────────┐
│ Cartesian   │────▶│ Domain      │
│  Topology   │     │ Decomposition│
└─────────────┘     └──────┬──────┘
                            │
                            ▼
                    ┌─────────────┐
                    │   Mesh2D   │
                    │  + Ghost   │
                    └──────┬──────┘
                           │
                           ▼
                    ┌─────────────┐     ┌─────────────┐
                    │   Initial   │────▶│   Time      │
                    │  Condition  │     │  Stepping   │
                    └─────────────┘     └──────┬──────┘
                                              │
                          ┌─────────────────────┼─────────────────────┐
                          │                     │                     │
                          ▼                     ▼                     ▼
                   ┌─────────────┐     ┌─────────────┐     ┌─────────────┐
                   │   Jacobi   │     │    SOR     │     │     CG     │
                   │   Solver   │     │   Solver   │     │   Solver   │
                   └──────┬──────┘     └──────┬──────┘     └──────┬──────┘
                          │                     │                     │
                          └─────────────────────┼─────────────────────┘
                                                │
                                                ▼
                                         ┌─────────────┐
                                         │   Output    │
                                         │  Exporter   │
                                         └─────────────┘
```

## Class Hierarchy

```
ISolver (Abstract)
├── JacobiSolver<Parallel>
├── SORSolver
└── ConjugateGradientSolver

CartesianTopology (Move-only)
├── Automatically computed dimensions
└── User-specified dimensions

Mesh2D (Copyable + Movable)
├── Serial version (no ghost cells)
└── Parallel version (with ghost cells)

DomainDecomposition (Move-only)
├── Automatic decomposition
├── Cartesian decomposition
└── Manual decomposition

CoordinateSystem (Copyable + Movable)
├── Uniform grid
└── Non-uniform grid
```

## Parallel Design

### Domain Decomposition Strategy

The solver uses 2D domain decomposition:

1. **Grid Partitioning**: Global domain divided into P subdomains
2. **Cartesian Topology**: Processes arranged in 2D grid
3. **Ghost Cells**: Each subdomain has overlapping boundary layers
4. **Communication**: Nearest neighbor exchange pattern

### Communication Patterns

```
┌────┬────┬────┐
│ P0 │ P1 │ P2 │  Exchange:
├────┼────┼────┤
│ P3 │ P4 │ P5 │  - North/South: MPI_Sendrecv
├────┼────┼────┤
│ P6 │ P7 │ P8 │  - East/West: MPI_Sendrecv
└────┴────┴────┘
```

### Parallel Algorithm Flow

1. **Initialize**: MPI, topology, domain decomposition
2. **Distribute**: Partition domain across processes
3. **Iterate**:
   - Exchange ghost cells (non-blocking)
   - Compute interior points (overlapping)
   - Compute residuals
   - Check convergence (global reduction)
4. **Gather**: Collect results on root process
5. **Output**: Write solution to file

### Load Balancing

Automatic dimension calculation ensures optimal load balance:
- Dimensions closest to square: minimize aspect ratio
- Even distribution: balanced workload
- Boundary processes: slightly less load (acceptable)

## Performance Considerations

### Memory Layout

- **Row-major storage**: Cache-friendly access patterns
- **Contiguous memory**: Efficient DMA transfers
- **Ghost cells**: Allocated for all processes (simpler code)

### Computation Overlap

- **Non-blocking communication**: Overlap exchange with computation
- **Ghost cell exchange**: Hidden latency
- **Asynchronous operations**: Better resource utilization

### Algorithm Selection

| Solver | Convergence | Parallelism | Best For |
|--------|-------------|-------------|----------|
| Jacobi | Slow | Excellent | Education, simple problems |
| SOR | Medium | Good (red-black) | Medium-sized problems |
| CG | Fast | Medium | Large-scale problems |

### Vectorization

- **SIMD-friendly**: Memory layout enables auto-vectorization
- **Loop structure**: Compiler can optimize inner loops
- **Stencil operations**: Good for CPU vector units

## Scalability Analysis

### Strong Scaling

Fixed problem size, increasing number of processes:

```
Speedup = T(1) / T(P)
Efficiency = Speedup / P

Ideal: Linear scaling
Real: Sub-linear due to communication overhead
```

### Weak Scaling

Fixed problem size per process, increasing total size:

```
Time per iteration ≈ constant
Communication overhead scales with P
```

### Communication Complexity

- **Per iteration**: O(√P) messages (Cartesian topology)
- **Message size**: O(N/P) (local domain size)
- **Latency**: O(log P) (for global reductions)

## Design Trade-offs

### Code Simplicity vs. Performance

**Choice**: Favor code simplicity with acceptable performance
- Clear, maintainable code
- Optimizations where beneficial
- Compiler can handle many optimizations

### Serial vs. Parallel Code

**Choice**: Template-based code sharing
- Single implementation for both
- Compile-time optimization
- Easy to maintain

### Memory vs. Computation

**Choice**: Trade memory for computation speed
- Ghost cells for fast neighbor access
- Pre-allocated work arrays
- Acceptable memory overhead

## Future Enhancements

### Potential Improvements

1. **Hybrid Parallelization**: MPI + OpenMP
2. **GPU Acceleration**: CUDA or OpenCL support
3. **Adaptive Mesh Refinement**: Dynamic grid resolution
4. **Multigrid Solver**: Faster convergence
5. **Load Balancing**: Dynamic domain redistribution
6. **Fault Tolerance**: Checkpoint/restart mechanism

### Extensibility Points

1. **New Solvers**: Implement ISolver interface
2. **New Time Schemes**: Extend time stepping layer
3. **New Mesh Types**: Implement mesh interface
4. **Custom Boundary Conditions**: Function-based BCs
5. **Custom Reductions**: Extend MPI operations

## References

- MPI Standard: https://www.mpi-forum.org/
- Modern C++ Design: Alexandrescu
- Parallel Programming for Multicore and Clusters: Gropp et al.
- Numerical Recipes: Press et al.
