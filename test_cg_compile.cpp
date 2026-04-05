/**
 * @file test_cg_compile.cpp
 * @brief Simple compilation test for CG solver
 */

#include <mpi.h>
#include <iostream>
#include "core/solver/conjugate_gradient_solver.hpp"
#include "utils/array2d.hpp"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        std::cout << "Testing CG solver compilation..." << std::endl;
    }

    // Create a simple problem
    int N = 30;
    double lambda = 0.1;

    utils::Array2D rhs(N, N, 1.0);
    utils::Array2D solution(N, N, 0.0);

    // Test standard CG
    {
        ConjugateGradientSolver solver(false, lambda);

        SolverParams params;
        params.tolerance = 1e-6;
        params.max_iterations = 100;
        params.verbose = (rank == 0);
        params.lambda = lambda;

        solver.solve(rhs, solution, params);

        if (rank == 0) {
            std::cout << "CG solver test passed!" << std::endl;
            std::cout << "  Iterations: " << solver.get_stats().iterations << std::endl;
            std::cout << "  Residual: " << solver.get_stats().residual << std::endl;
        }
    }

    // Test PCG
    solution.fill(0.0);
    {
        ConjugateGradientSolver solver(true, lambda);

        SolverParams params;
        params.tolerance = 1e-6;
        params.max_iterations = 100;
        params.verbose = (rank == 0);
        params.lambda = lambda;

        solver.solve(rhs, solution, params);

        if (rank == 0) {
            std::cout << "PCG solver test passed!" << std::endl;
            std::cout << "  Iterations: " << solver.get_stats().iterations << std::endl;
            std::cout << "  Residual: " << solver.get_stats().residual << std::endl;
        }
    }

    MPI_Finalize();
    return 0;
}
