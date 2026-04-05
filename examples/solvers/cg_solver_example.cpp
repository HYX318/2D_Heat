/**
 * @file cg_solver_example.cpp
 * @brief Example demonstrating Conjugate Gradient solver usage
 */

#include <mpi.h>
#include <iostream>
#include <iomanip>
#include "core/solver/conjugate_gradient_solver.hpp"
#include "utils/array2d.hpp"
#include <cmath>

int main(int argc, char** argv) {
    // Initialize MPI
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        std::cout << "=== Conjugate Gradient Solver Example ===" << std::endl;
        std::cout << "Number of processes: " << size << std::endl;
        std::cout << std::endl;
    }

    // Problem parameters
    int N = 50;              // Grid size (N x N)
    double lambda = 0.1;     // Diffusion coefficient dt/h^2

    // Create right-hand side vector
    utils::Array2D rhs(N, N, 0.0);
    utils::Array2D solution(N, N, 0.0);

    // Initialize RHS with a simple pattern (sinusoidal)
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            double x = static_cast<double>(i) / (N - 1);
            double y = static_cast<double>(j) / (N - 1);
            rhs(j, i) = std::sin(M_PI * x) * std::sin(M_PI * y);
        }
    }

    // Create solver parameters
    SolverParams params;
    params.tolerance = 1e-8;
    params.max_iterations = 500;
    params.verbose = (rank == 0);  // Only rank 0 prints
    params.lambda = lambda;

    // Run standard CG solver
    if (rank == 0) {
        std::cout << "--- Standard CG Solver ---" << std::endl;
        std::cout << "Grid size: " << N << " x " << N << std::endl;
        std::cout << "Lambda: " << lambda << std::endl;
        std::cout << std::endl;
    }

    {
        ConjugateGradientSolver solver(false, lambda);
        solver.solve(rhs, solution, params);

        const auto& stats = solver.get_stats();

        if (rank == 0) {
            std::cout << "CG Results:" << std::endl;
            std::cout << "  Iterations: " << stats.iterations << std::endl;
            std::cout << "  Residual: " << std::scientific << stats.residual << std::endl;
            std::cout << "  Solve time: " << std::fixed << stats.solve_time << " s" << std::endl;
            std::cout << "  Time/iter: " << stats.solve_time / stats.iterations << " s" << std::endl;
            std::cout << std::endl;
        }
    }

    // Run Preconditioned CG solver
    if (rank == 0) {
        std::cout << "--- Preconditioned CG Solver (Jacobi) ---" << std::endl;
    }

    solution.fill(0.0);  // Reset solution

    {
        ConjugateGradientSolver solver(true, lambda);
        solver.solve(rhs, solution, params);

        const auto& stats = solver.get_stats();

        if (rank == 0) {
            std::cout << "PCG Results:" << std::endl;
            std::cout << "  Iterations: " << stats.iterations << std::endl;
            std::cout << "  Residual: " << std::scientific << stats.residual << std::endl;
            std::cout << "  Solve time: " << std::fixed << stats.solve_time << " s" << std::endl;
            std::cout << "  Time/iter: " << stats.solve_time / stats.iterations << " s" << std::endl;
            std::cout << std::endl;
        }
    }

    // Compare with different lambda values
    if (rank == 0) {
        std::cout << "--- Performance Comparison ---" << std::endl;
        std::cout << std::setw(10) << "Lambda"
                  << std::setw(12) << "CG Iters"
                  << std::setw(12) << "PCG Iters"
                  << std::setw(15) << "CG Time (s)"
                  << std::setw(15) << "PCG Time (s)"
                  << std::endl;
        std::cout << std::string(65, '-') << std::endl;
    }

    std::vector<double> lambdas = {0.01, 0.05, 0.1, 0.2, 0.5};

    for (double lam : lambdas) {
        params.lambda = lam;
        params.verbose = false;

        // CG
        ConjugateGradientSolver cg_solver(false, lam);
        utils::Array2D sol_cg(N, N, 0.0);
        cg_solver.solve(rhs, sol_cg, params);
        const auto& cg_stats = cg_solver.get_stats();

        // PCG
        ConjugateGradientSolver pcg_solver(true, lam);
        utils::Array2D sol_pcg(N, N, 0.0);
        pcg_solver.solve(rhs, sol_pcg, params);
        const auto& pcg_stats = pcg_solver.get_stats();

        if (rank == 0) {
            std::cout << std::fixed << std::setprecision(2)
                      << std::setw(10) << lam
                      << std::setw(12) << cg_stats.iterations
                      << std::setw(12) << pcg_stats.iterations
                      << std::setprecision(4)
                      << std::setw(15) << cg_stats.solve_time
                      << std::setw(15) << pcg_stats.solve_time
                      << std::endl;
        }
    }

    if (rank == 0) {
        std::cout << std::endl;
        std::cout << "--- Summary ---" << std::endl;
        std::cout << "CG iterations ~ O(√κ) where κ is the condition number" << std::endl;
        std::cout << "For Poisson equation: κ = O(N²), so CG needs O(N) iterations" << std::endl;
        std::cout << "PCG (Jacobi preconditioner) typically reduces iterations by 2-3x" << std::endl;
    }

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
