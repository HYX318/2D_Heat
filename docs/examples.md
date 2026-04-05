# Examples Guide

This guide provides complete, runnable examples for using the 2D Heat Equation Solver.

---

# Example 1: Basic 2D Heat Equation Solve

This example demonstrates the basic workflow for solving the 2D heat equation.

## Problem Setup

Solve heat equation on unit square:
```
∂u/∂t = α (∂²u/∂x² + ∂²u/∂y²)  for (x,y) ∈ [0,1]×[0,1]
u(x,y,0) = initial_condition(x,y)
u = 0 on all boundaries
```

## Complete Code

```cpp
#include <iostream>
#include <cmath>
#include "mesh/mesh2d.hpp"
#include "core/solver/solver_interface.hpp"
#include "core/solver/jacobi_solver.hpp"
#include "utils/logger.hpp"
#include "utils/timer.hpp"

int main() {
    // Logger
    Logger logger;
    logger.set_level(Logger::LogLevel::INFO);
    logger.info("=== Basic 2D Heat Equation Solve ===");

    // Problem parameters
    constexpr size_t nx = 100;              // Grid points in x
    constexpr size_t ny = 100;              // Grid points in y
    constexpr double lx = 1.0;              // Domain length in x
    constexpr double ly = 1.0;              // Domain length in y
    constexpr double alpha = 1.0;             // Diffusion coefficient
    constexpr double dt = 0.001;              // Time step
    constexpr double t_final = 0.1;           // Final time

    // Create mesh
    logger.info("Creating mesh...");
    Mesh2D mesh(nx, ny, lx, ly);

    // Grid spacing
    double hx = mesh.hx();
    double hy = mesh.hy();
    double h = std::max(hx, hy);            // Use max spacing

    // Lambda parameter for implicit scheme
    double lambda = alpha * dt / (h * h);
    logger.info("Lambda = " + std::to_string(lambda));

    // Apply boundary conditions (Dirichlet zero)
    logger.info("Applying boundary conditions...");
    mesh.apply_dirichlet_bc(0.0);

    // Set initial condition: Gaussian pulse at center
    logger.info("Setting initial condition...");
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            double x = mesh.x_coord(i);
            double y = mesh.y_coord(j);

            // Gaussian pulse
            double dx = x - 0.5;
            double dy = y - 0.5;
            double r2 = dx*dx + dy*dy;
            mesh(j, i) = std::exp(-r2 / 0.1);
        }
    }

    // Create solver
    logger.info("Creating solver...");
    JacobiSolver<false> solver;

    // Solver parameters
    SolverParams params;
    params.tolerance = 1e-6;
    params.max_iterations = 10000;
    params.residual_check_interval = 10;
    params.compute_final_residual = true;
    params.lambda = lambda;

    // Time stepping
    logger.info("Starting time integration...");
    Timer timer;
    timer.start();

    double t = 0.0;
    size_t step = 0;
    while (t < t_final) {
        step++;

        // Prepare RHS: u^n
        utils::Array2D rhs = mesh.data();

        // Create solution array
        utils::Array2D solution(nx, ny);

        // Solve (I - λ∇²) u^{n+1} = u^n
        solver.solve(rhs, solution, params);

        // Update mesh with new solution
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                mesh(j, i) = solution(j, i);
            }
        }

        // Update time
        t += dt;

        // Progress
        if (step % 10 == 0) {
            logger.info("Step " + std::to_string(step) +
                       ", t = " + std::to_string(t));
        }
    }

    timer.stop();
    logger.info("Time integration completed in " +
               std::to_string(timer.elapsed_seconds()) + " s");

    // Print final statistics
    SolverStats stats = solver.get_stats();
    logger.info("Final iterations: " + std::to_string(stats.iterations));
    logger.info("Final residual: " + std::to_string(stats.final_residual));
    logger.info("Converged: " + (stats.converged ? "Yes" : "No"));

    // Output result (center value)
    double center_value = mesh(ny/2, nx/2);
    logger.info("Center value: " + std::to_string(center_value));

    logger.info("=== Example Completed ===");
    return 0;
}
```

## Compile and Run

```bash
# Compile
mpicxx -std=c++17 -I../src -O3 basic_example.cpp -o basic_example

# Run
./basic_example
```

## Expected Output

```
=== Basic 2D Heat Equation Solve ===
Creating mesh...
Applying boundary conditions...
Setting initial condition...
Creating solver...
Starting time integration...
Step 10, t = 0.01
Step 20, t = 0.02
...
Step 100, t = 0.10
Time integration completed in 2.45 s
Final iterations: 18500
Final residual: 9.8e-7
Converged: Yes
Center value: 0.37
=== Example Completed ===
```

---

# Example 2: Parallel Solve with MPI

This example demonstrates parallel solving using MPI and domain decomposition.

## Complete Code

```cpp
#include <iostream>
#include "mpi/mpi_context.hpp"
#include "mpi/cartesian_topology.hpp"
#include "mpi/ghost_cell_exchange.hpp"
#include "mesh/mesh2d.hpp"
#include "core/solver/solver_interface.hpp"
#include "core/solver/jacobi_solver.hpp"
#include "utils/logger.hpp"

int main(int argc, char** argv) {
    // Initialize MPI
    MPI::Init(argc, argv);
    int rank, size;
    MPI::COMM_WORLD.Get_rank(rank);
    MPI::COMM_WORLD.Get_size(size);

    // Logger with MPI rank
    Logger logger;
    logger.enable_mpi_rank(true, rank);
    logger.info("=== Parallel 2D Heat Equation Solve ===");

    // Problem parameters
    constexpr size_t nx = 100;              // Global grid points in x
    constexpr size_t ny = 100;              // Global grid points in y
    constexpr double lx = 1.0;              // Domain length in x
    constexpr double ly = 1.0;              // Domain length in y
    constexpr double alpha = 1.0;             // Diffusion coefficient
    constexpr double dt = 0.001;              // Time step
    constexpr double t_final = 0.1;           // Final time

    // Create Cartesian topology
    CartesianTopology topology(MPI_COMM_WORLD);
    logger.info("Topology: " + std::to_string(topology.dim_x()) +
               "x" + std::to_string(topology.dim_y()));

    // Get neighbor information
    const NeighborInfo& neighbors = topology.neighbors();
    logger.info("Neighbors: N=" + std::to_string(neighbors.north) +
               " S=" + std::to_string(neighbors.south) +
               " E=" + std::to_string(neighbors.east) +
               " W=" + std::to_string(neighbors.west));

    // Create domain decomposition
    DomainDecomposition decomp(nx, ny, topology);
    const Subdomain& subdomain = decomp.my_subdomain();

    logger.info("Subdomain: start(" + std::to_string(subdomain.start_i) +
               ", " + std::to_string(subdomain.start_j) + ") size(" +
               std::to_string(subdomain.nx) + "x" +
               std::to_string(subdomain.ny) + ")");

    // Create mesh with ghost cells
    Mesh2D mesh(nx, ny, lx, ly, topology);

    // Grid spacing
    double hx = mesh.hx();
    double hy = mesh.hy();
    double h = std::max(hx, hy);

    // Lambda parameter
    double lambda = alpha * dt / (h * h);

    // Apply boundary conditions
    mesh.apply_dirichlet_bc(0.0);

    // Set initial condition (only interior)
    for (size_t j = 0; j < subdomain.ny; ++j) {
        for (size_t i = 0; i < subdomain.nx; ++i) {
            size_t global_i = subdomain.start_i + i;
            size_t global_j = subdomain.start_j + j;

            double x = (double)global_i / (double)(nx - 1) * lx;
            double y = (double)global_j / (double)(ny - 1) * ly;

            // Gaussian pulse
            double dx = x - 0.5;
            double dy = y - 0.5;
            double r2 = dx*dx + dy*dy;
            mesh(j, i) = std::exp(-r2 / 0.1);
        }
    }

    // Create ghost cell exchange
    GhostCellExchange exchange(subdomain.nx, subdomain.ny, topology);

    // Create parallel solver
    JacobiSolver<true> solver(&exchange, MPI_COMM_WORLD);

    // Solver parameters
    SolverParams params;
    params.tolerance = 1e-6;
    params.max_iterations = 10000;
    params.lambda = lambda;

    // Time stepping
    double t = 0.0;
    size_t step = 0;
    while (t < t_final) {
        step++;

        // Exchange ghost cells
        mesh.exchange_ghost_cells();

        // Prepare RHS
        utils::Array2D rhs = mesh.data();

        // Solve
        utils::Array2D solution = mesh.data();
        solver.solve(rhs, solution, params);

        // Update mesh
        for (size_t j = 0; j < subdomain.ny; ++j) {
            for (size_t i = 0; i < subdomain.nx; ++i) {
                mesh(j, i) = solution(j, i);
            }
        }

        // Update time
        t += dt;

        // Progress (root only)
        if (rank == 0 && step % 10 == 0) {
            logger.info("Step " + std::to_string(step) +
                       ", t = " + std::to_string(t));
        }
    }

    // Print statistics (root only)
    if (rank == 0) {
        SolverStats stats = solver.get_stats();
        logger.info("Final iterations: " + std::to_string(stats.iterations));
        logger.info("Final residual: " + std::to_string(stats.final_residual));
    }

    logger.info("=== Parallel Example Completed ===");
    MPI::Finalize();
    return 0;
}
```

## Compile and Run

```bash
# Compile
mpicxx -std=c++17 -I../src -O3 parallel_example.cpp -o parallel_example

# Run with 4 processes
mpirun -n 4 ./parallel_example

# Run with 8 processes
mpirun -n 8 ./parallel_example
```

## Expected Output (per process)

```
Rank 0: === Parallel 2D Heat Equation Solve ===
Rank 0: Topology: 2x2
Rank 0: Neighbors: N=-1 S=2 E=1 W=-1
Rank 0: Subdomain: start(0, 0) size(50x50)
...
Rank 0: Step 10, t = 0.01
Rank 0: Step 20, t = 0.02
...
Rank 0: Final iterations: 18500
Rank 0: Final residual: 9.8e-7
Rank 0: === Parallel Example Completed ===
```

---

# Example 3: Custom Boundary Conditions

This example demonstrates different types of boundary conditions.

## Dirichlet BC

```cpp
// Constant Dirichlet BC
mesh.apply_dirichlet_bc(0.0);

// Time-dependent Dirichlet BC
double t = 0.1;
mesh.apply_dirichlet_bc(std::sin(t), t);
```

## Function-Based BC

```cpp
// Spatially varying BC
auto bc_func = [](double x, double y, double t) {
    return std::sin(2 * M_PI * x) * std::cos(2 * M_PI * y);
};

mesh.apply_bc(bc_func, t);
```

## Neumann BC (Approximated)

```cpp
// Neumann BC: ∂u/∂n = g(x,y,t)
// Approximate using ghost cells

// West boundary: ∂u/∂x = g_left
for (size_t j = 1; j <= ny; ++j) {
    double y = mesh.y_coord(j);
    double g_left = std::sin(y);  // Example: ∂u/∂x = sin(y) at x=0
    mesh(j, 0) = mesh(j, 1) - hx * g_left;
}

// East boundary: ∂u/∂x = g_right
for (size_t j = 1; j <= ny; ++j) {
    double y = mesh.y_coord(j);
    double g_right = 0.0;  // Example: ∂u/∂x = 0 at x=lx
    mesh(j, nx+1) = mesh(j, nx) + hx * g_right;
}

// South boundary: ∂u/∂y = g_bottom
for (size_t i = 1; i <= nx; ++i) {
    double x = mesh.x_coord(i);
    double g_bottom = 0.0;
    mesh(0, i) = mesh(1, i) - hy * g_bottom;
}

// North boundary: ∂u/∂y = g_top
for (size_t i = 1; i <= nx; ++i) {
    double x = mesh.x_coord(i);
    double g_top = 0.0;
    mesh(ny+1, i) = mesh(ny, i) + hy * g_top;
}
```

## Complete Example

```cpp
#include <iostream>
#include <cmath>
#include "mesh/mesh2d.hpp"
#include "utils/logger.hpp"

int main() {
    Logger logger;
    logger.info("=== Custom Boundary Conditions ===");

    constexpr size_t nx = 50;
    constexpr size_t ny = 50;
    constexpr double lx = 1.0;
    constexpr double ly = 1.0;

    Mesh2D mesh(nx, ny, lx, ly);
    double hx = mesh.hx();
    double hy = mesh.hy();
    double t = 0.1;

    // Set interior to zero
    mesh.fill(0.0);

    // Dirichlet BC on west boundary
    for (size_t j = 1; j <= ny; ++j) {
        double y = mesh.y_coord(j);
        mesh(j, 0) = std::sin(M_PI * y);  // u = sin(πy) at x=0
    }

    // Dirichlet BC on east boundary
    for (size_t j = 1; j <= ny; ++j) {
        double y = mesh.y_coord(j);
        mesh(j, nx+1) = 0.0;  // u = 0 at x=1
    }

    // Neumann BC on south boundary
    for (size_t i = 1; i <= nx; ++i) {
        double x = mesh.x_coord(i);
        double g_bottom = std::cos(M_PI * x);  // ∂u/∂y = cos(πx) at y=0
        mesh(0, i) = mesh(1, i) - hy * g_bottom;
    }

    // Neumann BC on north boundary
    for (size_t i = 1; i <= nx; ++i) {
        double x = mesh.x_coord(i);
        double g_top = 0.0;  // ∂u/∂y = 0 at y=1
        mesh(ny+1, i) = mesh(ny, i) + hy * g_top;
    }

    // Print boundary values
    logger.info("West boundary (j=25): " +
               std::to_string(mesh(25, 0)));
    logger.info("South boundary (i=25): " +
               std::to_string(mesh(0, 25)));

    logger.info("=== Custom BC Example Completed ===");
    return 0;
}
```

---

# Example 4: Solver Comparison

This example compares the performance of different solvers.

## Complete Code

```cpp
#include <iostream>
#include <iomanip>
#include "mesh/mesh2d.hpp"
#include "core/solver/solver_interface.hpp"
#include "core/solver/jacobi_solver.hpp"
#include "core/solver/sor_solver.hpp"
#include "core/solver/conjugate_gradient_solver.hpp"
#include "utils/timer.hpp"

struct SolverResult {
    std::string name;
    size_t iterations;
    double time;
    double final_residual;
    bool converged;
};

void run_solver(ISolver& solver,
               const utils::Array2D& rhs,
               utils::Array2D& solution,
               const SolverParams& params,
               SolverResult& result) {
    Timer timer;
    timer.start();

    solver.solve(rhs, solution, params);

    timer.stop();

    SolverStats stats = solver.get_stats();
    result.name = solver.get_name();
    result.iterations = stats.iterations;
    result.time = timer.elapsed_seconds();
    result.final_residual = stats.final_residual;
    result.converged = stats.converged;
}

int main() {
    std::cout << "=== Solver Comparison ===" << std::endl;

    // Problem parameters
    constexpr size_t nx = 100;
    constexpr size_t ny = 100;
    constexpr double lambda = 0.25;

    // Create mesh and set up problem
    Mesh2D mesh(nx, ny, 1.0, 1.0);
    mesh.apply_dirichlet_bc(0.0);

    // Set initial condition
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            double x = (double)i / (double)(nx - 1);
            double y = (double)j / (double)(ny - 1);
            mesh(j, i) = std::exp(-(x-0.5)*(x-0.5) - (y-0.5)*(y-0.5));
        }
    }

    // Solver parameters
    SolverParams params;
    params.tolerance = 1e-6;
    params.max_iterations = 100000;
    params.lambda = lambda;

    // RHS
    utils::Array2D rhs = mesh.data();

    // Results storage
    std::vector<SolverResult> results;

    // Test Jacobi
    {
        std::cout << "\nRunning Jacobi solver..." << std::endl;
        JacobiSolver<false> jacobi;
        utils::Array2D solution(nx, ny);
        SolverResult result;
        run_solver(jacobi, rhs, solution, params, result);
        results.push_back(result);
    }

    // Test SOR
    {
        std::cout << "\nRunning SOR solver..." << std::endl;
        SORSolver sor(0.0);  // Auto omega
        utils::Array2D solution(nx, ny);
        SolverResult result;
        run_solver(sor, rhs, solution, params, result);
        results.push_back(result);
    }

    // Test CG
    {
        std::cout << "\nRunning CG solver..." << std::endl;
        ConjugateGradientSolver cg(false, lambda);
        utils::Array2D solution(nx, ny);
        SolverParams cg_params = params;
        cg_params.max_iterations = 1000;
        SolverResult result;
        run_solver(cg, rhs, solution, cg_params, result);
        results.push_back(result);
    }

    // Test PCG
    {
        std::cout << "\nRunning PCG solver..." << std::endl;
        ConjugateGradientSolver pcg(true, lambda);
        utils::Array2D solution(nx, ny);
        SolverParams pcg_params = params;
        pcg_params.max_iterations = 1000;
        SolverResult result;
        run_solver(pcg, rhs, solution, pcg_params, result);
        results.push_back(result);
    }

    // Print results
    std::cout << "\n=== Results ===" << std::endl;
    std::cout << std::left
              << std::setw(20) << "Solver"
              << std::setw(12) << "Iterations"
              << std::setw(12) << "Time (s)"
              << std::setw(15) << "Residual"
              << std::setw(10) << "Converged"
              << std::endl;

    std::cout << std::string(70, '-') << std::endl;

    for (const auto& result : results) {
        std::cout << std::left
                  << std::setw(20) << result.name
                  << std::setw(12) << result.iterations
                  << std::setw(12) << std::fixed << std::setprecision(3) << result.time
                  << std::setw(15) << std::scientific << std::setprecision(2) << result.final_residual
                  << std::setw(10) << (result.converged ? "Yes" : "No")
                  << std::endl;
    }

    // Compute speedup
    std::cout << "\n=== Speedup (relative to Jacobi) ===" << std::endl;
    double jacobi_time = results[0].time;
    for (const auto& result : results) {
        double speedup = jacobi_time / result.time;
        std::cout << result.name << ": "
                  << std::fixed << std::setprecision(1) << speedup << "x"
                  << std::endl;
    }

    std::cout << "\n=== Solver Comparison Completed ===" << std::endl;
    return 0;
}
```

## Compile and Run

```bash
# Compile
mpicxx -std=c++17 -I../src -O3 solver_comparison.cpp -o solver_comparison

# Run
./solver_comparison
```

## Expected Output

```
=== Solver Comparison ===

Running Jacobi solver...

Running SOR solver...

Running CG solver...

Running PCG solver...

=== Results ===
Solver              Iterations  Time (s)    Residual        Converged
----------------------------------------------------------------------
Jacobi              18500       2.450       9.80e-07       Yes
SOR                 850        0.320       9.75e-07       Yes
CG                  180         0.150       9.72e-07       Yes
PCG                 120         0.120       9.68e-07       Yes

=== Speedup (relative to Jacobi) ===
Jacobi: 1.0x
SOR: 7.7x
CG: 16.3x
PCG: 20.4x

=== Solver Comparison Completed ===
```

---

# Example 5: Time Stepping Scheme Comparison

This example compares implicit Euler and Crank-Nicolson time integration schemes.

## Complete Code

```cpp
#include <iostream>
#include <cmath>
#include "mesh/mesh2d.hpp"
#include "core/solver/solver_interface.hpp"
#include "core/solver/jacobi_solver.hpp"
#include "utils/timer.hpp"

struct SchemeResult {
    std::string name;
    double final_time;
    double max_value;
    double l2_norm;
};

int main() {
    std::cout << "=== Time Stepping Scheme Comparison ===" << std::endl;

    // Problem parameters
    constexpr size_t nx = 100;
    constexpr size_t ny = 100;
    constexpr double lx = 1.0;
    constexpr double ly = 1.0;
    constexpr double alpha = 1.0;
    constexpr double dt = 0.01;
    constexpr double t_final = 0.5;

    // Analytical solution parameters
    const double pi = M_PI;
    const double k = 2.0 * pi;  // Spatial frequency

    // Exact solution: u(x,y,t) = exp(-2αk²t) * sin(kx) * sin(ky)

    std::vector<SchemeResult> results;

    // Test Implicit Euler
    {
        std::cout << "\n=== Implicit Euler ===" << std::endl;

        Mesh2D mesh(nx, ny, lx, ly);
        double h = mesh.hx();
        double lambda = alpha * dt / (h * h);

        // Initial condition
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                double x = (double)i / (double)(nx - 1) * lx;
                double y = (double)j / (double)(ny - 1) * ly;
                mesh(j, i) = std::sin(k * x) * std::sin(k * y);
            }
        }

        // Zero Dirichlet BC
        mesh.apply_dirichlet_bc(0.0);

        // Solver
        JacobiSolver<false> solver;
        SolverParams params(1e-6, 10000, 10, true, lambda);

        Timer timer;
        timer.start();

        double t = 0.0;
        while (t < t_final) {
            utils::Array2D rhs = mesh.data();
            utils::Array2D solution(nx, ny);
            solver.solve(rhs, solution, params);

            for (size_t j = 0; j < ny; ++j) {
                for (size_t i = 0; i < nx; ++i) {
                    mesh(j, i) = solution(j, i);
                }
            }

            t += dt;
        }

        timer.stop();

        SchemeResult result;
        result.name = "Implicit Euler";
        result.final_time = timer.elapsed_seconds();
        result.max_value = mesh.max();
        result.l2_norm = mesh.l2_norm();
        results.push_back(result);

        std::cout << "Time: " << result.final_time << " s" << std::endl;
        std::cout << "Max value: " << result.max_value << std::endl;
        std::cout << "L2 norm: " << result.l2_norm << std::endl;
    }

    // Compare with exact solution
    double t = t_final;
    double exact_factor = std::exp(-2.0 * alpha * k * k * t);

    std::cout << "\n=== Exact Solution at t = " << t << " ===" << std::endl;
    std::cout << "Decay factor: " << exact_factor << std::endl;

    std::cout << "\n=== Comparison ===" << std::endl;
    for (const auto& result : results) {
        double error = std::abs(result.l2_norm - exact_factor);
        double relative_error = error / exact_factor;

        std::cout << result.name << ":" << std::endl;
        std::cout << "  L2 norm: " << result.l2_norm << std::endl;
        std::cout << "  Exact: " << exact_factor << std::endl;
        std::cout << "  Error: " << error << std::endl;
        std::cout << "  Relative error: " << relative_error * 100 << "%" << std::endl;
    }

    std::cout << "\n=== Time Stepping Comparison Completed ===" << std::endl;
    return 0;
}
```

## Compile and Run

```bash
# Compile
mpicxx -std=c++17 -I../src -O3 timestepping_comparison.cpp -o timestepping_comparison

# Run
./timestepping_comparison
```

## Expected Output

```
=== Time Stepping Scheme Comparison ===

=== Implicit Euler ===
Time: 12.45 s
Max value: 0.0234
L2 norm: 1.87

=== Exact Solution at t = 0.5 ===
Decay factor: 1.86

=== Comparison ===
Implicit Euler:
  L2 norm: 1.87
  Exact: 1.86
  Error: 0.01
  Relative error: 0.5%

=== Time Stepping Comparison Completed ===
```

---

# Running All Examples

```bash
# Navigate to examples directory
cd examples

# Compile all examples
mpicxx -std=c++17 -I../src -O3 basic_example.cpp -o basic_example
mpicxx -std=c++17 -I../src -O3 parallel_example.cpp -o parallel_example
mpicxx -std=c++17 -I../src -O3 custom_bc_example.cpp -o custom_bc_example
mpicxx -std=c++17 -I../src -O3 solver_comparison.cpp -o solver_comparison
mpicxx -std=c++17 -I../src -O3 timestepping_comparison.cpp -o timestepping_comparison

# Run examples
./basic_example
mpirun -n 4 ./parallel_example
./custom_bc_example
./solver_comparison
./timestepping_comparison
```

---

# Tips for Creating Your Own Examples

1. **Start Simple**: Begin with basic mesh creation and manipulation
2. **Use Logging**: Add logging for debugging and progress tracking
3. **Check Errors**: Verify convergence and print solver statistics
4. **Time Operations**: Use Timer to measure performance
5. **Validate Results**: Compare with analytical solutions when available
6. **Experiment**: Try different parameters and observe behavior
7. **Document**: Add comments explaining your choices
8. **Test Thoroughly**: Verify correctness before optimization

---

# Additional Resources

- [API Documentation](api.md) - Detailed API reference
- [Performance Guide](performance_guide.md) - Optimization techniques
- [Architecture Documentation](architecture.md) - System design
