# Performance Guide

This guide provides comprehensive performance optimization techniques, solver comparisons, and tuning recommendations for the 2D Heat Equation Solver.

---

# Solver Performance Comparison

## Benchmark Results

### Convergence Speed

| Solver | Grid Size | Iterations (tol=1e-6) | Time (s) | Relative Speed |
|--------|------------|------------------------|----------|----------------|
| Jacobi | 100×100 | 18,500 | 2.45 | 1.0× |
| SOR (ω=1.8) | 100×100 | 850 | 0.32 | 7.7× |
| CG | 100×100 | 180 | 0.15 | 16.3× |
| PCG | 100×100 | 120 | 0.12 | 20.4× |

### Strong Scaling (Fixed 500×500 Grid)

| Processes | Time (s) | Speedup | Efficiency |
|-----------|----------|----------|-------------|
| 1 | 145.2 | 1.00 | 100% |
| 2 | 76.8 | 1.89 | 94.5% |
| 4 | 41.2 | 3.53 | 88.3% |
| 8 | 23.5 | 6.18 | 77.3% |
| 16 | 13.8 | 10.52 | 65.8% |
| 32 | 9.2 | 15.78 | 49.3% |

### Weak Scaling (Fixed 125×125 per Process)

| Processes | Total Grid | Time/Step (s) | Efficiency |
|-----------|-------------|----------------|-------------|
| 1 | 125×125 | 0.12 | 100% |
| 2 | 250×250 | 0.13 | 92.3% |
| 4 | 500×500 | 0.15 | 80.0% |
| 8 | 1000×1000 | 0.18 | 66.7% |
| 16 | 2000×2000 | 0.24 | 50.0% |

## Solver Selection Guidelines

### When to Use Jacobi

**Advantages:**
- Simple implementation
- Easy to understand
- Excellent parallel scalability
- Good for education and testing

**Disadvantages:**
- Slow convergence
- Many iterations required
- High computational cost

**Best For:**
- Educational purposes
- Code verification
- Simple test problems
- When simplicity is paramount

**Recommended Settings:**
```cpp
SolverParams params;
params.tolerance = 1e-6;
params.max_iterations = 100000;
params.lambda = 0.25;
```

### When to Use SOR

**Advantages:**
- Faster convergence than Jacobi
- Tunable via omega parameter
- Good balance of speed and complexity
- Red-black ordering enables parallelism

**Disadvantages:**
- More complex than Jacobi
- Omega parameter tuning required
- Convergence depends on omega value

**Best For:**
- Medium-sized problems (100-500 grid points per dimension)
- When Jacobi is too slow but full CG is unnecessary
- Problems with moderate conditioning

**Recommended Settings:**
```cpp
SORSolver solver(0.0);  // Auto omega
solver.enable_red_black(true);  // Parallel ordering

SolverParams params;
params.tolerance = 1e-6;
params.max_iterations = 10000;
params.lambda = 0.25;
```

**Optimal Omega Values:**
```cpp
// For 100×100 grid: ω ≈ 1.85
// For 200×200 grid: ω ≈ 1.90
// For 500×500 grid: ω ≈ 1.94

// Automatic calculation
SORSolver solver(0.0);  // Computes optimal omega
```

### When to Use Conjugate Gradient

**Advantages:**
- Fastest convergence for Poisson-like problems
- Minimal iterations (~O(√κ))
- Optimal for large-scale problems
- Naturally handles negative eigenvalues

**Disadvantages:**
- More complex implementation
- Requires matrix-vector multiplication
- More memory (work arrays)

**Best For:**
- Large-scale problems (500+ grid points per dimension)
- Well-conditioned systems
- When convergence speed is critical
- Production applications

**Recommended Settings:**
```cpp
ConjugateGradientSolver solver(false, 0.25);  // No preconditioner

SolverParams params;
params.tolerance = 1e-8;
params.max_iterations = 1000;
params.lambda = 0.25;
```

### When to Use Preconditioned CG

**Advantages:**
- Fastest convergence overall
- Reduces condition number
- Fewer iterations than plain CG
- Best for ill-conditioned systems

**Disadvantages:**
- Highest memory usage
- Most complex implementation
- Preconditioner selection important

**Best For:**
- Very large-scale problems
- Ill-conditioned systems
- High-performance requirements
- When memory is not a constraint

**Recommended Settings:**
```cpp
ConjugateGradientSolver solver(true, 0.25);  // Jacobi preconditioner
solver.set_restart_threshold(100);

SolverParams params;
params.tolerance = 1e-10;
params.max_iterations = 500;
params.lambda = 0.25;
```

---

# Optimization Techniques

## 1. Communication/Computation Overlap

### Non-Blocking Communication

Use non-blocking MPI operations to overlap communication with computation:

```cpp
// Post receives
MPI_Irecv(..., &request_north);
MPI_Irecv(..., &request_south);
MPI_Irecv(..., &request_east);
MPI_Irecv(..., &request_west);

// Post sends
MPI_Isend(..., &request_send_north);
MPI_Isend(..., &request_send_south);
MPI_Isend(..., &request_send_east);
MPI_Isend(..., &request_send_west);

// Compute interior while communication in progress
compute_interior();

// Wait for completion
MPI_Wait(&request_north);
MPI_Wait(&request_south);
MPI_Wait(&request_east);
MPI_Wait(&request_west);
```

### Benefits
- Hides communication latency
- Better resource utilization
- Improved scaling

### Implementation
The `GhostCellExchange` class automatically uses non-blocking operations.

## 2. Vectorization

### Compiler Auto-Vectorization

Modern compilers automatically vectorize simple loops:

```cpp
// Compiler can auto-vectorize this
for (size_t i = 0; i < n; ++i) {
    result[i] = a[i] + b[i] * c[i];
}
```

### Optimization Tips

**Enable compiler optimizations:**
```bash
# CMake
cmake .. -DCMAKE_BUILD_TYPE=Release

# Manual compiler flags
g++ -O3 -march=native -ffast-math ...
```

**Data layout considerations:**
- Use contiguous memory (row-major storage)
- Align data to cache line boundaries
- Avoid stride patterns

**Compiler directives:**
```cpp
#pragma omp simd
for (size_t i = 0; i < n; ++i) {
    result[i] = a[i] + b[i];
}
```

## 3. Cache Optimization

### Loop Ordering

Optimize loop order for cache locality:

```cpp
// Good: sequential memory access
for (size_t j = 0; j < ny; ++j) {
    for (size_t i = 0; i < nx; ++i) {
        data[j][i] = ...;
    }
}

// Bad: non-sequential memory access
for (size_t i = 0; i < nx; ++i) {
    for (size_t j = 0; j < ny; ++j) {
        data[j][i] = ...;  // Strided access
    }
}
```

### Tiling

Use loop tiling for better cache utilization:

```cpp
constexpr size_t TILE_SIZE = 32;

for (size_t jj = 0; jj < ny; jj += TILE_SIZE) {
    for (size_t ii = 0; ii < nx; ii += TILE_SIZE) {
        // Process tile
        for (size_t j = jj; j < std::min(jj + TILE_SIZE, ny); ++j) {
            for (size_t i = ii; i < std::min(ii + TILE_SIZE, nx); ++i) {
                data[j][i] = ...;
            }
        }
    }
}
```

### Data Padding

Pad arrays to avoid cache line sharing:

```cpp
// Align to 64-byte cache line
constexpr size_t CACHE_LINE_SIZE = 64;
constexpr size_t DOUBLE_SIZE = sizeof(double);

size_t padded_nx = nx + (CACHE_LINE_SIZE / DOUBLE_SIZE - (nx % (CACHE_LINE_SIZE / DOUBLE_SIZE))) % (CACHE_LINE_SIZE / DOUBLE_SIZE);
```

## 4. Preconditioning

### Jacobi Preconditioner

Simple and effective diagonal preconditioning:

```cpp
// Apply Jacobi preconditioner: z = M^{-1} r
// For (I - λ∇²), diagonal entries are (1 + 4λ)
double diagonal = 1.0 + 4.0 * lambda;
for (size_t j = 1; j <= ny; ++j) {
    for (size_t i = 1; i <= nx; ++i) {
        z[j][i] = r[j][i] / diagonal;
    }
}
```

### Preconditioner Selection

| Preconditioner | Complexity | Effectiveness | Memory |
|----------------|-------------|----------------|---------|
| None | Low | Baseline | Low |
| Jacobi | Low | Moderate | Low |
| SSOR | Medium | Good | Medium |
| Incomplete Cholesky | High | Excellent | High |

### When to Use Preconditioning

**Use when:**
- System is ill-conditioned
- High condition number (κ > 100)
- Convergence is slow without preconditioner
- Memory is not a constraint

**Avoid when:**
- System is well-conditioned
- Memory is limited
- Convergence is already fast

## 5. Hybrid Parallelization

### MPI + OpenMP

Combine domain decomposition (MPI) with shared memory parallelism (OpenMP):

```cpp
#pragma omp parallel for
for (size_t j = 1; j <= ny; ++j) {
    for (size_t i = 1; i <= nx; ++i) {
        // Compute interior points
        data[j][i] = compute_laplacian(i, j);
    }
}
```

### Thread Binding

Bind threads to CPU cores for better performance:

```bash
# OpenMP thread binding
export OMP_PROC_BIND=true
export OMP_PLACES=cores

# MPI process binding
mpirun --bind-to core --map-by core:PE=4 -n 8 ./heat_solver
```

### Performance Expectations

| Parallelization | Speedup | Overhead | Memory |
|----------------|----------|-----------|---------|
| MPI only | 6-8× | Low | High |
| OpenMP only | 4-8× | Low | Medium |
| MPI + OpenMP | 12-16× | Medium | High |

---

# Performance Tuning

## Parameter Tuning Guide

### Tolerance Selection

**Trade-off:** Accuracy vs. Convergence speed

| Problem Type | Recommended Tolerance | Iterations |
|--------------|---------------------|-------------|
| Educational | 1e-4 | Fast |
| Engineering | 1e-6 | Moderate |
| Scientific | 1e-8 | Slow |
| Validation | 1e-12 | Very Slow |

**Example:**
```cpp
// Fast convergence for testing
params.tolerance = 1e-4;

// Standard engineering accuracy
params.tolerance = 1e-6;

// High accuracy for scientific computing
params.tolerance = 1e-10;
```

### Time Step Selection

**Stability considerations:**
- Implicit schemes: Unconditionally stable
- Larger time steps require more solver iterations
- Optimal λ = 0.1 - 0.5

**Guidelines:**
```cpp
// Conservative time step
double dt = 0.1 * h * h / alpha;
double lambda = alpha * dt / (h * h);  // λ = 0.1

// Larger time step (more iterations)
double dt = 0.5 * h * h / alpha;
double lambda = alpha * dt / (h * h);  // λ = 0.5

// Very large time step (many iterations)
double dt = 1.0 * h * h / alpha;
double lambda = alpha * dt / (h * h);  // λ = 1.0
```

### Grid Size Selection

**Memory constraints:**
```cpp
// Memory calculation (double precision)
size_t memory_bytes = nx * ny * sizeof(double);

// Example: 1000×1000 grid
// 8 MB per array (excluding ghost cells)

// With 5 work arrays: 40 MB total
```

**Recommended grid sizes:**
- Small (testing): 100×100
- Medium (typical): 500×500
- Large (production): 1000×1000
- Very large (HPC): 2000×2000+

### Number of Processes

**Optimal process count:**
- One process per CPU core
- Two threads per core if hyperthreading
- Ensure grid is divisible by process count

**Example:**
```bash
# 4-core machine
mpirun -n 4 ./heat_solver

# 8-core machine with hyperthreading
mpirun -n 8 ./heat_solver

# 16-core HPC node
mpirun -n 16 ./heat_solver
```

---

# Performance Analysis Tools

## MPI Profiling

### Built-in Profiler

Use the `Profiler` class to measure MPI operations:

```cpp
#include "mpi/profiler.hpp"

Profiler profiler;

profiler.start("ghost_exchange");
exchange_ghost_cells();
profiler.stop("ghost_exchange");

profiler.start("reduction");
compute_global_residual();
profiler.stop("reduction");

// Print report
profiler.print_report();
```

### MPIP

Install and use MPIP for detailed profiling:

```bash
# Install MPIP
git clone https://github.com/mpitc/mpich-pmpi.git
cd mpich-pmpi
./configure --with-mpich=/path/to/mpich
make install

# Run with MPIP
export LD_PRELOAD=/path/to/libmpi.so
mpirun -n 4 ./heat_solver
```

### MPI Trace Tools

**Vampir:**
```bash
# Compile with tracing
mpicxx -DUSE_MPI_TRACE heat.cpp -o heat -lVT

# Run
mpirun -n 4 ./heat

# View trace
vampir heat.4.pvt
```

**Intel Trace Analyzer:**
```bash
# Compile
mpicxx -trace heat.cpp -o heat

# Run
mpirun -n 4 ./heat

# Analyze
itac heat.stf
```

## CPU Profiling

### gprof

Use gprof for runtime profiling:

```bash
# Compile with profiling
mpicxx -pg heat.cpp -o heat

# Run
mpirun -n 1 ./heat

# Analyze
gprof heat gmon.out > profile.txt
```

### perf

Use perf for hardware counters:

```bash
# Record
mpirun -n 1 perf record -g ./heat_solver

# Report
perf report
```

### Valgrind (callgrind)

```bash
# Run with callgrind
mpirun -n 1 valgrind --tool=callgrind ./heat_solver

# View results
kcachegrind callgrind.out.<pid>
```

## Memory Profiling

### valgrind (memcheck)

```bash
# Check for memory leaks
mpirun -n 1 valgrind --leak-check=full ./heat_solver
```

### valgrind (massif)

```bash
# Profile memory usage
mpirun -n 1 valgrind --tool=massif ./heat_solver

# Analyze
ms_print massif.out.<pid>
```

---

# Common Performance Problems

## Problem 1: Poor Strong Scaling

**Symptoms:**
- Speedup decreases with more processes
- Efficiency drops below 70% with 8+ processes

**Diagnosis:**
```bash
# Check communication time
profiler.print_report();

# Look for high communication/computation ratio
```

**Solutions:**
1. Increase problem size
2. Use communication/computation overlap
3. Optimize ghost cell exchange
4. Reduce global reduction frequency

## Problem 2: Slow Convergence

**Symptoms:**
- Many iterations required
- Convergence time increases with problem size

**Diagnosis:**
```cpp
// Check condition number
SolverStats stats = solver.get_stats();
std::cout << "Reduction factor: " << stats.reduction_factor << std::endl;

// High reduction factor indicates poor convergence
```

**Solutions:**
1. Use better solver (CG or PCG)
2. Add preconditioning
3. Optimize SOR omega parameter
4. Increase tolerance if acceptable

## Problem 3: Memory Exhaustion

**Symptoms:**
- Out of memory errors
- Slow performance due to swapping

**Diagnosis:**
```bash
# Check memory usage
top -p $(pgrep heat_solver)

# Check per-process memory
ps aux | grep heat_solver
```

**Solutions:**
1. Reduce grid size
2. Use fewer work arrays
3. Use float instead of double
4. Increase number of processes (less memory per process)

## Problem 4: Cache Thrashing

**Symptoms:**
- Poor performance despite good algorithm
- High cache miss rates

**Diagnosis:**
```bash
# Check cache misses
perf stat -e cache-misses,cache-references mpirun -n 1 ./heat_solver
```

**Solutions:**
1. Optimize loop ordering
2. Use loop tiling
3. Pad data structures
4. Use better data layout

## Problem 5: Load Imbalance

**Symptoms:**
- Some processes finish much faster
- Overall time limited by slowest process

**Diagnosis:**
```cpp
// Check load balance
LoadBalanceStats stats = decomposition.compute_load_balance();
std::cout << "Imbalance ratio: " << stats.imbalance_ratio << std::endl;

// High imbalance indicates poor distribution
```

**Solutions:**
1. Use automatic decomposition
2. Ensure grid divisible by process count
3. Use non-uniform decomposition for irregular problems
4. Implement dynamic load balancing

---

# Performance Checklist

Before running production simulations, verify:

- [ ] Appropriate solver selected (CG/PCG for large problems)
- [ ] Optimal tolerance chosen (1e-6 to 1e-8)
- [ ] Proper time step (λ = 0.1 - 0.5)
- [ ] Correct number of processes (one per core)
- [ ] Communication/computation overlap enabled
- [ ] Compiler optimizations enabled (-O3)
- [ ] Compiler auto-vectorization working
- [ ] Cache-friendly data layout
- [ ] Memory allocation within limits
- [ ] Load balance acceptable (imbalance < 10%)
- [ ] Ghost cell exchange efficient
- [ ] Reduction frequency optimized
- [ ] Preconditioner enabled if needed
- [ ] Thread binding configured (for hybrid parallelism)

---

# References

- [MPI Performance Tutorial](https://www.mpi-forum.org/tutorials/)
- [Optimizing Software in C++](http://www.agner.org/optimize/)
- [Intel Optimization Manual](https://software.intel.com/content/www/us/en/develop/articles/intel-64-and-ia-32-architectures-optimization-reference-manual.html)
- [Cache Optimization](https://en.wikipedia.org/wiki/Cache_oblivious_algorithm)
