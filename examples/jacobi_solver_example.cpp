/**
 * @file jacobi_solver_example.cpp
 * @brief Example demonstrating the use of the refactored Jacobi solver
 *
 * This example shows:
 * 1. How to create and configure the solver
 * 2. How to set up a simple test problem
 * 3. How to solve the system and retrieve results
 * 4. How to use parallel solver with MPI
 */

#include <iostream>
#include <iomanip>
#include <mpi.h>
#include <cmath>

#include "src/core/solver/solver_interface.hpp"
#include "src/core/solver/jacobi_solver.hpp"
#include "src/utils/array2d.hpp"
#include "src/mpi/ghost_cell_exchange.hpp"
#include "src/mpi/cartesian_topology.hpp"

using namespace utils;

/**
 * Create a simple test problem with known analytical solution
 *
 * Problem: -Δu = f on unit square with u = 0 on boundaries
 * Exact solution: u(x,y) = sin(πx) × sin(πy)
 */
void create_test_problem(int Nx, int Ny, double lambda,
                       Array2D& rhs, Array2D& exact) {
    double pi = 3.14159265358979323846;

    for (int j = 0; j < Ny + 2; ++j) {
        for (int i = 0; i < Nx + 2; ++i) {
            double x = static_cast<double>(i) / (Nx + 1);
            double y = static_cast<double>(j) / (Ny + 1);

            // Exact solution
            exact(j, i) = std::sin(pi * x) * std::sin(pi * y);

            // For interior points, compute RHS
            if (i >= 1 && i <= Nx && j >= 1 && j <= Ny) {
                double laplacian = exact(j+1, i) - 2.0*exact(j, i) + exact(j-1, i) +
                                exact(j, i+1) - 2.0*exact(j, i) + exact(j, i-1);
                rhs(j, i) = exact(j, i) - lambda * laplacian;
            } else {
                // Boundary conditions
                rhs(j, i) = exact(j, i);
            }
        }
    }
}

/**
 * Compute maximum error between computed and exact solution
 */
double compute_max_error(const Array2D& computed, const Array2D& exact,
                       int Nx, int Ny) {
    double max_error = 0.0;

    for (int j = 1; j <= Ny; ++j) {
        for (int i = 1; i <= Nx; ++i) {
            double error = std::abs(computed(j, i) - exact(j, i));
            max_error = std::max(max_error, error);
        }
    }

    return max_error;
}

/**
 * Serial solver example
 */
void serial_solver_example() {
    std::cout << "\n=== Serial Jacobi Solver Example ===\n" << std::endl;

    // Problem parameters
    int Nx = 50;  // Interior columns
    int Ny = 50;  // Interior rows
    double lambda = 0.25;

    std::cout << "Problem size: " << Nx << " x " << Ny << std::endl;
    std::cout << "Lambda: " << lambda << std::endl;

    // Create arrays (including ghost cells)
    Array2D rhs(Ny + 2, Nx + 2);
    Array2D exact(Ny + 2, Nx + 2);
    Array2D solution(Ny + 2, Nx + 2);

    // Create test problem
    create_test_problem(Nx, Ny, lambda, rhs, exact);

    // Create solver
    JacobiSolver<false> solver;

    // Configure solver parameters
    SolverParams params;
    params.tolerance = 1e-8;
    params.max_iterations = 10000;
    params.lambda = lambda;

    // Solve
    std::cout << "\nSolving..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    solver.solve(rhs, solution, params);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    double solve_time = elapsed.count();

    // Get statistics
    SolverStats stats = solver.get_stats();

    // Compute error
    double max_error = compute_max_error(solution, exact, Nx, Ny);

    // Print results
    std::cout << "\nResults:" << std::endl;
    std::cout << "  Converged: " << (stats.converged ? "Yes" : "No") << std::endl;
    std::cout << "  Iterations: " << stats.iterations << std::endl;
    std::cout << "  Final residual: " << stats.final_residual << std::endl;
    std::cout << "  Solve time: " << solve_time << " seconds" << std::endl;
    std::cout << "  Time per iteration: " << solve_time / stats.iterations << " seconds" << std::endl;
    std::cout << "  Max error: " << max_error << std::endl;
}

/**
 * Parallel solver example
 */
void parallel_solver_example() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        std::cout << "\n=== Parallel Jacobi Solver Example ===\n" << std::endl;
        std::cout << "Number of processes: " << size << std::endl;
    }

    // Problem parameters
    int global_Nx = 50;  // Total interior columns
    int global_Ny = 50;  // Total interior rows
    double lambda = 0.25;

    // Create Cartesian topology
    CartesianTopology topology(MPI_COMM_WORLD);

    // Each process handles a subdomain
    int dim_x = topology.dim_x();
    int dim_y = topology.dim_y();

    int proc_x = topology.coord_x();
    int proc_y = topology.coord_y();

    int nx_per_proc = global_Nx / dim_x;
    int ny_per_proc = global_Ny / dim_y;

    if (rank == 0) {
        std::cout << "Topology: " << dim_x << " x " << dim_y << std::endl;
        std::cout << "Problem size: " << global_Nx << " x " << global_Ny << std::endl;
        std::cout << "Subdomain per process: " << nx_per_proc << " x " << ny_per_proc << std::endl;
        std::cout << "Lambda: " << lambda << std::endl;
    }

    // Create local arrays (including ghost cells)
    Array2D local_rhs(ny_per_proc + 2, nx_per_proc + 2);
    Array2D local_exact(ny_per_proc + 2, nx_per_proc + 2);
    Array2D local_solution(ny_per_proc + 2, nx_per_proc + 2);

    // Fill local problem with test data
    // (In a real application, this would be computed from global problem)
    local_rhs.fill(0.0);
    local_exact.fill(0.0);

    // Add some simple pattern for demonstration
    for (int j = 1; j <= ny_per_proc; ++j) {
        for (int i = 1; i <= nx_per_proc; ++i) {
            int global_i = proc_x * nx_per_proc + i;
            int global_j = proc_y * ny_per_proc + j;

            double x = static_cast<double>(global_i) / (global_Nx + 1);
            double y = static_cast<double>(global_j) / (global_Ny + 1);
            double pi = 3.14159265358979323846;

            local_exact(j, i) = std::sin(pi * x) * std::sin(pi * y);
            local_rhs(j, i) = local_exact(j, i) * (1.0 + 4.0 * lambda);  // Simplified RHS
        }
    }

    // Create ghost cell exchange
    GhostCellExchange ghost_exchange(nx_per_proc, ny_per_proc, topology);

    // Create parallel solver
    JacobiSolver<true> solver(&ghost_exchange, MPI_COMM_WORLD);

    // Configure solver parameters
    SolverParams params;
    params.tolerance = 1e-8;
    params.max_iterations = 10000;
    params.lambda = lambda;

    // Solve
    MPI_Barrier(MPI_COMM_WORLD);
    auto start = std::chrono::high_resolution_clock::now();

    solver.solve(local_rhs, local_solution, params);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    double solve_time = elapsed.count();

    // Get statistics (same on all processes)
    SolverStats stats = solver.get_stats();

    // Print results (only rank 0)
    if (rank == 0) {
        std::cout << "\nSolving..." << std::endl;
        std::cout << "\nResults:" << std::endl;
        std::cout << "  Converged: " << (c ? "Yes" : "No") << std::endl;
        std::cout << "  Iterations: " << stats.iterations << std::endl;
        std::cout << "  Final residual: " << stats.final_residual << std::endl;
        std::cout << "  Solve time: " << solve_time << " seconds" << std::endl;
        std::cout << "  Time per iteration: " << solve_time / stats.iterations << " seconds" << std::endl;
    }
}

int main(int argc, char** argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::cout << std::fixed << std::setprecision(8);

    // Run serial example (only on rank 0)
    if (rank == 0) {
        serial_solver_example();
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Run parallel example (on all ranks)
    parallel_solver_example();

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
