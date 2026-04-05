# API Documentation

This document provides detailed API documentation for all public interfaces in the 2D Heat Equation Solver.

---

# MPI Layer

## MPIContext

### Description
RAII wrapper for MPI initialization and finalization. Provides safe, exception-safe MPI management.

### Header
```cpp
#include "mpi/mpi_context.hpp"
```

### Public Methods

#### Constructor
```cpp
MPIContext(int& argc, char** argv);
```
**Parameters:**
- `argc`: Argument count from main()
- `argv`: Argument vector from main()

**Throws:** `std::runtime_error` if MPI initialization fails

**Example:**
```cpp
int main(int argc, char** argv) {
    MPIContext mpi(argc, argv);
    // MPI is now initialized
    std::cout << "Rank: " << mpi.rank() << std::endl;
    return 0;
} // MPI automatically finalized
```

#### Get Process Rank
```cpp
int rank() const;
```
**Returns:** Process rank (0 to size-1)

#### Get Number of Processes
```cpp
int size() const;
```
**Returns:** Total number of processes

#### Check if Root
```cpp
bool is_root() const;
```
**Returns:** `true` if rank == 0, `false` otherwise

#### Barrier Synchronization
```cpp
void barrier() const;
```
Blocks until all processes reach this point.

**Throws:** `std::runtime_error` if the barrier call fails

#### Abort Execution
```cpp
void abort(int error_code, const std::string& message) const;
```
**Parameters:**
- `error_code`: Error code to return to the environment
- `message`: Error message to display

Calls `MPI_Abort` to terminate all MPI processes.

---

## CartesianTopology

### Description
Manages MPI Cartesian topology for domain decomposition. Provides automatic dimension calculation, coordinate mapping, and neighbor identification.

### Header
```cpp
#include "mpi/cartesian_topology.hpp"
```

### Public Methods

#### Construct with Automatic Dimensions
```cpp
explicit CartesianTopology(MPI_Comm comm = MPI_COMM_WORLD);
```
Automatically calculates optimal dimensions (closest to square) based on number of processes.

**Throws:** `std::runtime_error` if topology creation fails

**Example:**
```cpp
CartesianTopology topology(MPI_COMM_WORLD);
std::cout << "Dimensions: " << topology.dim_x() << "x" << topology.dim_y() << std::endl;
```

#### Construct with Specified Dimensions
```cpp
CartesianTopology(MPI_Comm comm, const std::vector<int>& dims);
```
**Parameters:**
- `comm`: MPI communicator
- `dims`: Dimensions [dim_x, dim_y]

**Throws:**
- `std::invalid_argument` if dims.size() != 2
- `std::invalid_argument` if product of dims != number of processes
- `std::runtime_error` if topology creation fails

#### Get Communicator
```cpp
MPI_Comm communicator() const;
```
**Returns:** MPI communicator with Cartesian topology

#### Get Process Coordinates
```cpp
const std::vector<int>& coords() const;
```
**Returns:** Vector [coord_x, coord_y]

#### Get Neighbors
```cpp
const NeighborInfo& neighbors() const;
```
**Returns:** `NeighborInfo` structure with neighbor ranks

**NeighborInfo structure:**
```cpp
struct NeighborInfo {
    int south;  // Neighbor in -Y direction
    int north;  // Neighbor in +Y direction
    int west;   // Neighbor in -X direction
    int east;   // Neighbor in +X direction
};
```

#### Check Boundary Status
```cpp
bool is_on_boundary(Direction dir, Shift shift) const;
```
**Parameters:**
- `dir`: Direction (X or Y)
- `shift`: Shift direction (Forward or Backward)

**Returns:** `true` if on boundary, `false` otherwise

---

## GhostCellExchange

### Description
Efficient ghost cell exchange between MPI processes. Supports non-blocking communication and overlapping with computation.

### Header
```cpp
#include "mpi/ghost_cell_exchange.hpp"
```

### Public Methods

#### Constructor
```cpp
GhostCellExchange(size_t nx, size_t ny, const CartesianTopology& topology);
```
**Parameters:**
- `nx`: Local number of interior points in x-direction
- `ny`: Local number of interior points in y-direction
- `topology`: Cartesian topology for neighbor information

#### Exchange Ghost Cells
```cpp
void exchange(utils::Array2D& array, MPI_Comm comm = MPI_COMM_WORLD);
```
**Parameters:**
- `array`: Array with ghost cells (includes boundary layers)
- `comm`: MPI communicator (default: MPI_COMM_WORLD)

**Throws:** `std::runtime_error` if communication fails

**Example:**
```cpp
GhostCellExchange exchange(nx, ny, topology);
utils::Array2D data(ny+2, nx+2);

// After updating interior
exchange.exchange(data, MPI_COMM_WORLD);

// Ghost cells now contain neighbor data
```

---

## ReductionOps

### Description
Type-safe wrapper for MPI reduction operations. Supports scalar and array reductions, scan operations, and broadcast.

### Header
```cpp
#include "mpi/reduction_ops.hpp"
```

### Public Methods

#### Constructor
```cpp
ReductionOps();  // Uses MPI_COMM_WORLD
explicit ReductionOps(MPI_Comm comm);  // Uses custom communicator
explicit ReductionOps(const MPIContext& ctx);  // Uses MPIContext's communicator
```

#### Reduce Scalar to Root
```cpp
template<typename T>
T reduce(T value, ReductionOp op, int root) const;
```
**Parameters:**
- `value`: Value to reduce from this process
- `op`: Reduction operation (SUM, MAX, MIN, etc.)
- `root`: Rank of root process

**Returns:** Reduced value on root process, undefined on others

**Example:**
```cpp
ReductionOps ops;
int my_sum = compute_local_sum();
int global_sum = ops.reduce(my_sum, ReductionOp::SUM, 0);
```

#### Allreduce (Broadcast to All)
```cpp
template<typename T>
T allreduce(T value, ReductionOp op) const;
```
**Returns:** Reduced value (same on all processes)

#### Scan (Inclusive Prefix Reduction)
```cpp
template<typename T>
T scan(T value, ReductionOp op) const;
```
**Returns:** Reduction on processes 0 through current rank

#### Exclusive Scan
```cpp
template<typename T>
T exscan(T value, ReductionOp op) const;
```
**Returns:** Reduction on processes 0 through (current rank - 1)

#### Find Maximum and Location
```cpp
template<typename T>
MaxLocResult<T> maxloc(T value) const;
```
**Returns:** `MaxLocResult` containing max value and rank where it was found

**Result structure:**
```cpp
template<typename T>
struct MaxLocResult {
    T value;  // Maximum value
    int rank; // Rank where max was found
};
```

#### Broadcast
```cpp
template<typename T>
void broadcast(T& value, int root) const;
```
Broadcasts value from root to all processes.

---

## Profiler

### Description
Performance profiling for MPI operations. Collects timing and statistics for communication patterns.

### Header
```cpp
#include "mpi/profiler.hpp"
```

### Public Methods

#### Start Timing
```cpp
void start(const std::string& operation);
```
**Parameters:**
- `operation`: Name of operation to profile

#### Stop Timing
```cpp
void stop(const std::string& operation);
```
**Parameters:**
- `operation`: Name of operation being profiled

#### Get Statistics
```cpp
ProfilingStats get_stats(const std::string& operation) const;
```
**Returns:** Statistics for specified operation

#### Print Report
```cpp
void print_report() const;
```
Prints profiling report to stdout.

---

# Mesh Layer

## Mesh2D

### Description
2D rectangular mesh with ghost cell support. Handles coordinate transformations, boundary conditions, and numerical operations.

### Header
```cpp
#include "mesh/mesh2d.hpp"
```

### Public Methods

#### Constructor (Serial)
```cpp
Mesh2D(size_t nx, size_t ny, double lx = 1.0, double ly = 1.0);
```
**Parameters:**
- `nx`: Number of interior points in x-direction
- `ny`: Number of interior points in y-direction
- `lx`: Physical length in x-direction
- `ly`: Physical length in y-direction

**Throws:** `std::invalid_argument` if nx or ny is 0

#### Constructor (Parallel)
```cpp
Mesh2D(size_t nx, size_t ny, double lx, double ly,
       const CartesianTopology& topology);
```
Creates mesh with ghost cells for MPI domain decomposition.

#### Element Access
```cpp
double& operator()(size_t i, size_t j);
const double& operator()(size_t i, size_t j) const;
```
**Parameters:**
- `i`: Row index (0 to ny-1 for interior)
- `j`: Column index (0 to nx-1 for interior)

**Returns:** Reference to element at (i, j)

**Example:**
```cpp
Mesh2D mesh(100, 100, 1.0, 1.0);
mesh(50, 50) = 1.0;  // Set center value
double value = mesh(50, 50);  // Get center value
```

#### Get Physical Coordinates
```cpp
double x_coord(size_t i) const;
double y_coord(size_t j) const;
std::pair<double, double> coord(size_t i, size_t j) const;
```
**Returns:** Physical coordinate(s) at grid point

#### Apply Boundary Conditions
```cpp
void apply_dirichlet_bc(double value = 0.0, double t = 0.0);
void apply_bc(std::function<double(double x, double y, double t)> bc_func, double t = 0.0);
```
**Parameters:**
- `value`: Constant boundary value
- `bc_func`: Function defining boundary condition f(x, y, t)
- `t`: Time parameter

**Example:**
```cpp
// Constant BC
mesh.apply_dirichlet_bc(0.0);

// Function-based BC
auto bc_func = [](double x, double y, double t) {
    return std::sin(x) * std::cos(y);
};
mesh.apply_bc(bc_func, 0.0);
```

#### Compute Laplacian
```cpp
void compute_laplacian(utils::Array2D& result) const;
```
**Parameters:**
- `result`: Output array for Laplacian

Computes Laplacian using 5-point stencil:
```
∇²u = (u_{i+1,j} - 2u_{i,j} + u_{i-1,j})/hx² +
      (u_{i,j+1} - 2u_{i,j} + u_{i,j-1})/hy²
```

#### Norm Calculations
```cpp
double l2_norm() const;        // Euclidean norm
double linfty_norm() const;    // Maximum absolute value
double max() const;            // Maximum value
double min() const;            // Minimum value
```

#### Exchange Ghost Cells (Parallel)
```cpp
void exchange_ghost_cells();
```
Exchanges ghost cell data with neighboring processes.

---

## DomainDecomposition

### Description
Manages domain decomposition for parallel computing. Divides global domain into subdomains for MPI processes.

### Header
```cpp
#include "mesh/domain_decomposition.hpp"
```

### Public Methods

#### Construct with Automatic Decomposition
```cpp
DomainDecomposition(size_t global_nx, size_t global_ny,
                   const CartesianTopology& topology);
```
**Parameters:**
- `global_nx`: Global domain size in X direction
- `global_ny`: Global domain size in Y direction
- `topology`: Cartesian topology

Automatically calculates optimal decomposition.

#### Get Current Subdomain
```cpp
const Subdomain& my_subdomain() const;
```
**Returns:** Information about current process's subdomain

**Subdomain structure:**
```cpp
struct Subdomain {
    size_t start_i;      // Global starting index in X
    size_t start_j;      // Global starting index in Y
    size_t nx;           // Local domain size in X
    size_t ny;           // Local domain size in Y
    size_t global_nx;    // Global domain size in X
    size_t global_ny;    // Global domain size in Y
};
```

#### Coordinate Transformation
```cpp
bool is_in_my_domain(size_t global_i, size_t global_j) const;
int find_owner_rank(size_t global_i, size_t global_j) const;
```
Checks if coordinates are in current domain and finds owning process.

#### Load Balance Analysis
```cpp
LoadBalanceStats compute_load_balance() const;
```
**Returns:** Statistics about load balancing

---

## CoordinateSystem

### Description
Handles coordinate system, grid generation, and coordinate transformations. Supports uniform and non-uniform grids.

### Header
```cpp
#include "mesh/coordinate_system.hpp"
```

### Public Methods

#### Constructor (Uniform Grid)
```cpp
CoordinateSystem(size_t nx, size_t ny,
                double x_min = 0.0, double x_max = 1.0,
                double y_min = 0.0, double y_max = 1.0);
```
**Parameters:**
- `nx`, `ny`: Number of points in X and Y directions
- `x_min`, `x_max`: X coordinate range
- `y_min`, `y_max`: Y coordinate range

#### Coordinate Calculation
```cpp
double x(size_t i) const;
double y(size_t j) const;
std::pair<double, double> coord(size_t i, size_t j) const;
```
**Returns:** Physical coordinate(s) at grid index

**Example:**
```cpp
CoordinateSystem coords(100, 100, 0.0, 1.0, 0.0, 1.0);
double x_center = coords.x(50);  // 0.5
double y_center = coords.y(50);  // 0.5
```

#### Index Lookup
```cpp
size_t i_at(double x) const;
size_t j_at(double y) const;
```
**Returns:** Nearest grid index for physical coordinate

#### Non-Uniform Grid Support
```cpp
void set_x_distribution(std::function<double(double)> dist);
void set_y_distribution(std::function<double(double)> dist);
```
**Parameters:**
- `dist`: Function mapping normalized [0,1] to physical coordinate

**Example:**
```cpp
CoordinateSystem coords(100, 100);

// Boundary layer stretching
coords.set_x_distribution([](double xi) {
    double beta = 1.2;
    return CoordinateSystem::stretching_function(xi, beta);
});
```

#### Grid Generation
```cpp
void generate_coordinate_mesh(utils::Array2D& x_mesh,
                             utils::Array2D& y_mesh) const;
```
Generates coordinate meshgrids for visualization.

---

# Solver Layer

## ISolver (Interface)

### Description
Abstract base class for iterative solvers. Defines common interface for all solver implementations.

### Header
```cpp
#include "core/solver/solver_interface.hpp"
```

### Public Methods

#### Solve Linear System
```cpp
virtual void solve(const utils::Array2D& rhs,
                   utils::Array2D& solution,
                   const SolverParams& params) = 0;
```
**Parameters:**
- `rhs`: Right-hand side array
- `solution`: Solution array (used as initial guess, overwritten)
- `params`: Solver parameters

#### Get Statistics
```cpp
virtual SolverStats get_stats() const = 0;
```
**Returns:** Statistics from last solve operation

**SolverStats structure:**
```cpp
struct SolverStats {
    bool converged;          // Converged flag
    size_t iterations;      // Number of iterations
    double final_residual;  // Final residual norm
    double initial_residual; // Initial residual norm
    double solve_time;      // Solve time (seconds)
    double reduction_factor; // Residual reduction factor
};
```

#### Get Solver Name
```cpp
virtual std::string get_name() const = 0;
```
**Returns:** Solver name string

#### Reset Solver State
```cpp
virtual void reset() = 0;
```
Clears cached data and statistics.

---

## JacobiSolver

### Description
Jacobi iterative solver for 2D heat equation. Template-based for serial/parallel execution.

### Header
```cpp
#include "core/solver/jacobi_solver.hpp"
```

### Public Methods

#### Constructor (Serial)
```cpp
explicit JacobiSolver();
```

#### Constructor (Parallel)
```cpp
JacobiSolver(GhostCellExchange* ghost_exchange, MPI_Comm comm);
```
**Parameters:**
- `ghost_exchange`: Pointer to ghost cell exchange manager
- `comm`: MPI communicator

#### Solve
```cpp
void solve(const utils::Array2D& rhs,
           utils::Array2D& solution,
           const SolverParams& params) override;
```

**Example:**
```cpp
// Serial
JacobiSolver<false> solver;
solver.solve(rhs, solution, params);

// Parallel
GhostCellExchange exchange(nx, ny, topology);
JacobiSolver<true> solver(&exchange, MPI_COMM_WORLD);
solver.solve(rhs, solution, params);
```

---

## SORSolver

### Description
Successive Over-Relaxation solver with automatic omega optimization. Supports red-black ordering for parallelism.

### Header
```cpp
#include "core/solver/sor_solver.hpp"
```

### Public Methods

#### Constructor
```cpp
explicit SORSolver(double omega = 0.0);
SORSolver(double omega, GhostCellExchange* ghost_exchange, MPI_Comm comm = MPI_COMM_WORLD);
```
**Parameters:**
- `omega`: Relaxation parameter (0.0 for automatic calculation)
- `ghost_exchange`: Ghost cell exchange for parallel version
- `comm`: MPI communicator

**Omega values:**
- `0.0`: Automatically compute optimal omega
- `0.0 < omega < 1.0`: Under-relaxation (stable, slower)
- `1.0`: Gauss-Seidel (standard)
- `1.0 < omega < 2.0`: Over-relaxation (faster)
- `omega >= 2.0`: Invalid (will diverge)

#### Set Relaxation Parameter
```cpp
void set_omega(double omega);
```

#### Enable Red-Black Ordering
```cpp
void enable_red_black(bool enable);
```
Red-black ordering improves parallelization by allowing independent updates.

**Example:**
```cpp
SORSolver solver(0.0);  // Auto omega
solver.enable_red_black(true);
solver.solve(rhs, solution, params);
```

---

## ConjugateGradientSolver

### Description
Conjugate Gradient solver with optional Jacobi preconditioning. Superior convergence for Poisson-like problems.

### Header
```cpp
#include "core/solver/conjugate_gradient_solver.hpp"
```

### Public Methods

#### Constructor (Serial)
```cpp
explicit ConjugateGradientSolver(bool use_preconditioner = false, double lambda = 0.0);
```

#### Constructor (Parallel)
```cpp
ConjugateGradientSolver(bool use_preconditioner, double lambda,
                       MPI_Comm comm,
                       const int neighbor_rank[4],
                       int nx, int ny);
```

#### Solve
```cpp
void solve(const utils::Array2D& rhs,
           utils::Array2D& solution,
           const SolverParams& params);
```

#### Set Restart Threshold
```cpp
void set_restart_threshold(int max_iterations_without_progress);
```
Sets threshold for numerical stability restart mechanism.

**Example:**
```cpp
ConjugateGradientSolver solver(true, 0.25);  // PCG with lambda
solver.set_restart_threshold(100);
solver.solve(rhs, solution, params);
```

---

# Utilities Layer

## Array2D

### Description
RAII 2D array wrapper with automatic memory management, bounds checking, and mathematical operations.

### Header
```cpp
#include "utils/array2d.hpp"
```

### Public Methods

#### Constructor
```cpp
Array2D(size_t rows, size_t cols);                    // Uninitialized
Array2D(size_t rows, size_t cols, double value);     // Fill with value
```

#### Element Access
```cpp
double& operator()(size_t i, size_t j);
const double& operator()(size_t i, size_t j) const;
```

#### Dimensions
```cpp
size_t rows() const;
size_t cols() const;
size_t size() const;  // rows * cols
```

#### Fill and Copy
```cpp
void fill(double value);
void copy_from(const Array2D& other);
```

#### Norms
```cpp
double l2_norm() const;      // sqrt(sum of squares)
double linfty_norm() const;  // Maximum absolute value
double max() const;
double min() const;
```

#### Arithmetic Operations
```cpp
Array2D& operator+=(const Array2D& other);
Array2D& operator-=(const Array2D& other);
Array2D& operator*=(double scalar);
```

#### Raw Data Access
```cpp
double* data();
const double* data() const;
```
Returns pointer to contiguous memory storage.

**Example:**
```cpp
utils::Array2D grid(100, 100, 0.0);
grid(50, 50) = 1.0;
double norm = grid.l2_norm();
grid *= 2.0;  // Scale by 2
```

---

## Logger

### Description
Thread-safe structured logging with multiple levels, timestamps, and MPI rank support.

### Header
```cpp
#include "utils/logger.hpp"
```

### Public Methods

#### Constructor
```cpp
explicit Logger(LogLevel min_level = LogLevel::INFO);
Logger(const std::string& filename, LogLevel min_level = LogLevel::INFO);
```

#### Logging Methods
```cpp
void debug(const std::string& message);
void info(const std::string& message);
void warning(const std::string& message);
void error(const std::string& message);
void fatal(const std::string& message);
```

#### Configuration
```cpp
void set_level(LogLevel level);
void enable_timestamps(bool enable);
void enable_mpi_rank(bool enable, int rank = -1);
```

**Example:**
```cpp
Logger logger;
logger.enable_mpi_rank(true, mpi.rank());
logger.info("Starting simulation");
logger.error("Failed to open file");
```

---

## Timer

### Description
High-precision timer for performance measurement.

### Header
```cpp
#include "utils/timer.hpp"
```

### Public Methods

#### Start/Stop
```cpp
void start();
void stop();
```

#### Elapsed Time
```cpp
double elapsed_seconds() const;
double elapsed_milliseconds() const;
double elapsed_microseconds() const;
```

#### Reset
```cpp
void reset();
```

**Example:**
```cpp
Timer timer;
timer.start();
// ... computation ...
timer.stop();
std::cout << "Time: " << timer.elapsed_milliseconds() << " ms" << std::endl;
```

---

# Common Structures

## SolverParams

### Description
Parameters for iterative solver configuration.

```cpp
struct SolverParams {
    double tolerance;                // Convergence tolerance
    size_t max_iterations;           // Maximum iterations
    size_t residual_check_interval;   // Residual check frequency
    bool compute_final_residual;      // Compute final residual
    double lambda;                  // Time step coefficient

    SolverParams(double tol, size_t max_iter, size_t check_interval = 10,
                bool compute_residual = true, double lam = 0.25);
};
```

## Subdomain

### Description
Information about a local subdomain in domain decomposition.

```cpp
struct Subdomain {
    size_t start_i;      // Global starting index in X
    size_t start_j;      // Global starting index in Y
    size_t nx;           // Local domain size in X
    size_t ny;           // Local domain size in Y
    size_t global_nx;    // Global domain size in X
    size_t global_ny;    // Global domain size in Y
};
```

## NeighborInfo

### Description
Neighbor process ranks in Cartesian topology.

```cpp
struct NeighborInfo {
    int south;  // Neighbor in -Y direction
    int north;  // Neighbor in +Y direction
    int west;   // Neighbor in -X direction
    int east;   // Neighbor in +X direction
};
```

---

# Usage Examples

## Complete Serial Solver Example

```cpp
#include "mesh/mesh2d.hpp"
#include "core/solver/solver_interface.hpp"
#include "core/solver/jacobi_solver.hpp"
#include "utils/array2d.hpp"

int main() {
    // Create mesh
    Mesh2D mesh(100, 100, 1.0, 1.0);

    // Apply boundary conditions
    mesh.apply_dirichlet_bc(0.0);

    // Set initial condition
    mesh(50, 50) = 1.0;

    // Create solver
    Jac JacobiSolver solver;
    SolverParams params(1e-6, 10000, 10, true, 0.25);

    // Solve
    utils::Array2D rhs = mesh.data();
    utils::Array2D solution(100, 100);
    solver.solve(rhs, solution, params);

    // Get statistics
    SolverStats stats = solver.get_stats();
    std::cout << "Iterations: " << stats.iterations << std::endl;
    std::cout << "Converged: " << (stats.converged ? "Yes" : "No") << std::endl;

    return 0;
}
```

## Complete Parallel Solver Example

```cpp
#include "mpi/mpi_context.hpp"
#include "mpi/cartesian_topology.hpp"
#include "mpi/ghost_cell_exchange.hpp"
#include "mesh/mesh2d.hpp"
#include "core/solver/jacobi_solver.hpp"

int main(int argc, char** argv) {
    // Initialize MPI
    MPIContext mpi(argc, argv);

    // Create Cartesian topology
    CartesianTopology topology(MPI_COMM_WORLD);

    // Create mesh with ghost cells
    Mesh2D mesh(100, 100, 1.0, 1.0, topology);

    // Apply boundary conditions
    mesh.apply_dirichlet_bc(0.0);

    // Create ghost cell exchange
    GhostCellExchange exchange(mesh.nx(), mesh.ny(), topology);

    // Create parallel solver
    JacobiSolver<true> solver(&exchange, MPI_COMM_WORLD);
    SolverParams params(1e-6, 10000, 10, true, 0.25);

    // Solve
    utils::Array2D rhs = mesh.data();
    utils::Array2D solution = mesh.data();
    solver.solve(rhs, solution, params);

    // Print result on root
    if (mpi.is_root()) {
        SolverStats stats = solver.get_stats();
        std::cout << "Iterations: " << stats.iterations << std::endl;
    }

    return 0;
}
```
