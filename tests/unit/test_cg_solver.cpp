/**
 * @file test_cg_solver.cpp
 * @brief Unit tests for Conjugate Gradient solver
 */

#include <gtest/gtest.h>
#include "core/solver/conjugate_gradient_solver.hpp"
#include "utils/array2d.hpp"
#include <cmath>
#include <chrono>
#include <iostream>

using utils::Array2D;

/**
 * @test ConstructionTest - Test CG solver construction
 */
TEST(CGTest, ConstructionTest) {
    // Test serial construction
    {
        ConjugateGradientSolver solver(false, 0.1);
        EXPECT_FALSE(solver.is_parallel());
        EXPECT_EQ(solver.get_rank(), 0);
        EXPECT_EQ(solver.get_size(), 1);
        EXPECT_EQ(solver.get_type(), SolverType::ConjugateGradient);
        EXPECT_EQ(solver.get_restarts(), 0);
    }

    // Test with preconditioner
    {
        ConjugateGradientSolver solver(true, 0.1);
        EXPECT_EQ(solver.get_type(), SolverType::PreconditionedCG);
    }

    // Test parallel construction (mock)
    {
        int neighbor_rank[4] = {0, 0, 0, 0}; // Dummy neighbors
        ConjugateGradientSolver solver(false, 0.1, MPI_COMM_WORLD, neighbor_rank, 10, 10);
        EXPECT_TRUE(solver.is_parallel());
    }
}

/**
 * @test ConvergenceTest - Test CG convergence behavior
 */
TEST(CGTest, ConvergenceTest) {
    // Create test problem: (I - λ∇²)x = b
    int N = 50;
    double lambda = 0.1;
    Array2D rhs(N, N, 0.0);
    Array2D solution(N, N, 0.0);

    // Create a simple right-hand side
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            rhs(j, i) = 1.0;
        }
    }

    // Create solver
    ConjugateGradientSolver solver(false, lambda);

    // Set solver parameters
    SolverParams params;
    params.tolerance = 1e-8;
    params.max_iterations = 500;
    params.verbose = false;
    params.lambda = lambda;

    // Solve
    solver.solve(rhs, solution, params);

    // Get statistics
    const auto& stats = solver.get_stats();

    // Check convergence
    EXPECT_LT(stats.residual, params.tolerance);
    EXPECT_GT(stats.iterations, 0);
    EXPECT_LE(stats.iterations, params.max_iterations);
    EXPECT_GT(stats.solve_time, 0.0);

    // Verify solution is positive (for positive RHS)
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            EXPECT_GT(solution(j, i), 0.0);
        }
    }
}

/**
 * @test PCGConvergenceTest - Test Preconditioned CG convergence
 */
TEST(CGTest, PCGConvergenceTest) {
    int N = 50;
    double lambda = 0.1;
    Array2D rhs(N, N, 0.0);
    Array2D solution(N, N, 0.0);

    // Create RHS with a pattern
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            rhs(j, i) = std::sin(M_PI * i / N) * std::sin(M_PI * j / N);
        }
    }

    // Create PCG solver
    ConjugateGradientSolver solver(true, lambda);

    SolverParams params;
    params.tolerance = 1e-8;
    params.max_iterations = 500;
    params.verbose = false;
    params.lambda = lambda;

    solver.solve(rhs, solution, params);
    const auto& stats = solver.get_stats();

    EXPECT_LT(stats.residual, params.tolerance);
    EXPECT_EQ(solver.get_type(), SolverType::PreconditionedCG);
}

/**
 * @test SpeedupTest - Compare CG vs Jacobi performance
 */
TEST(CGTest, SpeedupTest) {
    int N = 100;
    double lambda = 0.05;

    // Create test problem
    Array2D rhs(N, N, 1.0);
    Array2D solution_cg(N, N, 0.0);
    Array2D solution_jacobi(N, N, 0.0);

    // CG solver
    ConjugateGradientSolver solver_cg(false, lambda);
    SolverParams params;
    params.tolerance = 1e-6;
    params.max_iterations = 1000;
    params.verbose = false;
    params.lambda = lambda;

    auto start_cg = std::chrono::high_resolution_clock::now();
    solver_cg.solve(rhs, solution_cg, params);
    auto end_cg = std::chrono::high_resolution_clock::now();
    double time_cg = std::chrono::duration<double>(end_cg - start_cg).count();

    // Jacobi solver (simplified implementation)
    auto start_jacobi = std::chrono::high_resolution_clock::now();

    Array2D x(N, N, 0.0);
    Array2D x_old(N, N, 0.0);
    double coeff = 1.0 / (1.0 + 4.0 * lambda);
    double residual = 1.0;
    int iter_jacobi = 0;

    while (residual > params.tolerance && iter_jacobi < params.max_iterations) {
        iter_jacobi++;
        residual = 0.0;

        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                double sum = 0.0;
                if (j > 0) sum += x_old(j - 1, i);
                if (j < N - 1) sum += x_old(j + 1, i);
                if (i > 0) sum += x_old(j, i - 1);
                if (i < N - 1) sum += x_old(j, i + 1);

                x(j, i) = coeff * (lambda * sum + rhs(j, i));
                residual = std::max(residual, std::abs(x(j, i) - x_old(j, i)));
            }
        }

        Array2D temp = std::move(x);
        x = std::move(x_old);
        x_old = std::move(temp);
    }

    auto end_jacobi = std::chrono::high_resolution_clock::now();
    double time_jacobi = std::chrono::duration<double>(end_jacobi - start_jacobi).count();

    // CG should be much faster
    std::cout << "\nCG iterations: " << solver_cg.get_stats().iterations
              << ", time: " << time_cg << " s" << std::endl;
    std::cout << "Jacobi iterations: " << iter_jacobi
              << ", time: " << time_jacobi << " s" << std::endl;

    // CG should converge in significantly fewer iterations
    EXPECT_LT(solver_cg.get_stats().iterations, iter_jacobi * 0.3);
}

/**
 * @test PreconditionerTest - Compare CG vs PCG
 */
TEST(CGTest, PreconditionerTest) {
    int N = 80;
    double lambda = 0.1;

    Array2D rhs(N, N, 0.0);
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            rhs(j, i) = std::cos(2.0 * M_PI * i / N) * std::sin(2.0 * M_PI * j / N);
        }
    }

    Array2D sol_cg(N, N, 0.0);
    Array2D sol_pcg(N, N, 0.0);

    // Standard CG
    ConjugateGradientSolver solver_cg(false, lambda);
    SolverParams params;
    params.tolerance = 1e-7;
    params.max_iterations = 1000;
    params.verbose = false;
    params.lambda = lambda;

    solver_cg.solve(rhs, sol_cg, params);

    // Preconditioned CG
    ConjugateGradientSolver solver_pcg(true, lambda);
    solver_pcg.solve(rhs, sol_pcg, params);

    // PCG should converge in fewer iterations
    EXPECT_LE(solver_pcg.get_stats().iterations, solver_cg.get_stats().iterations);

    // Both should converge to similar solution
    double diff = 0.0;
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            diff = std::max(diff, std::abs(sol_cg(j, i) - sol_pcg(j, i)));
        }
    }
    EXPECT_LT(diff, 1e-5);
}

/**
 * @test RestartTest - Test restart mechanism
 */
TEST(CGTest, RestartTest) {
    int N = 60;
    double lambda = 0.2;  // Larger lambda for more challenging problem

    Array2D rhs(N, N, 0.0);
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            rhs(j, i) = std::sin(M_PI * i / N) * std::cos(M_PI * j / N);
        }
    }

    Array2D solution(N, N, 0.0);

    // Create solver with low restart threshold
    ConjugateGradientSolver solver(false, lambda);
    solver.set_restart_threshold(10);

    SolverParams params;
    params.tolerance = 1e-6;
    params.max_iterations = 1000;
    params.verbose = false;
    params.lambda = lambda;

    solver.solve(rhs, solution, params);
    const auto& stats = solver.get_stats();

    // Should converge
    EXPECT_LT(stats.residual, params.tolerance);

    // Check if restarts occurred (may or may not depending on problem)
    std::cout << "\nRestarts: " << stats.restarts << std::endl;
}

/**
 * @test ParallelTest - Test parallel dot product correctness
 */
TEST(CGTest, ParallelTest) {
    int N = 40;
    double lambda = 0.1;

    Array2D rhs(N, N, 1.0);
    Array2D solution(N, N, 0.0);

    // Create parallel solver (even if running with 1 process)
    int neighbor_rank[4] = {-1, -1, -1, -1};

    ConjugateGradientSolver solver(false, lambda, MPI_COMM_WORLD, neighbor_rank, N, N);

    SolverParams params;
    params.tolerance = 1e-8;
    params.max_iterations = 500;
    params.verbose = false;
    params.lambda = lambda;

    solver.solve(rhs, solution, params);

    // Should converge
    EXPECT_LT(solver.get_stats().residual, params.tolerance);

    // Verify solution is positive
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            EXPECT_GT(solution(j, i), 0.0);
        }
    }
}

/**
 * @test ConditionNumberTest - Test behavior with different condition numbers
 */
TEST(CGTest, ConditionNumberTest) {
    // Different lambda values give different condition numbers
    std::vector<double> lambdas = {0.01, 0.1, 0.5, 1.0};

    std::cout << "\nCondition number test:" << std::endl;
    std::cout << "Lambda\tIterations\tResidual" << std::endl;

    for (double lambda : lambdas) {
        int N = 50;
        Array2D rhs(N, N, 1.0);
        Array2D solution(N, N, 0.0);

        ConjugateGradientSolver solver(false, lambda);
        SolverParams params;
        params.tolerance = 1e-6;
        params.max_iterations = 1000;
        params.verbose = false;
        params.lambda = lambda;

        solver.solve(rhs, solution, params);

        const auto& stats = solver.get_stats();

        std::cout << lambda << "\t" << stats.iterations
                  << "\t" << std::scientific << stats.residual << std::endl;

        // Should converge
        EXPECT_LT(stats.residual, params.tolerance);
    }
}

/**
 * @test ExceptionTest - Test error handling
 */
TEST(CGTest, ExceptionTest) {
    int N = 50;
    double lambda = 0.1;

    Array2D rhs(N, N, 1.0);
    Array2D solution(N, N, 0.0);

    ConjugateGradientSolver solver(false, lambda);
    SolverParams params;
    params.tolerance = 1e-8;
    params.max_iterations = 500;
    params.verbose = false;
    params.lambda = lambda;

    // Mismatched dimensions should throw
    Array2D wrong_size(N + 1, N, 0.0);
    EXPECT_THROW(solver.solve(rhs, wrong_size, params), std::invalid_argument);
}

/**
 * @test TrivialSolutionTest - Test zero RHS
 */
TEST(CGTest, TrivialSolutionTest) {
    int N = 30;
    double lambda = 0.1;

    Array2D rhs(N, N, 0.0);  // Zero RHS
    Array2D solution(N, N, 0.0);

    ConjugateGradientSolver solver(false, lambda);
    SolverParams params;
    params.tolerance = 1e-8;
    params.max_iterations = 500;
    params.verbose = false;
    params.lambda = lambda;

    solver.solve(rhs, solution, params);
    const auto& stats = solver.get_stats();

    // Should converge immediately
    EXPECT_LT(stats.residual, params.tolerance);
    EXPECT_EQ(stats.iterations, 0);  // Should be 0 for trivial solution
}

/**
 * @test ResetTest - Test reset functionality
 */
TEST(CGTest, ResetTest) {
    ConjugateGradientSolver solver(false, 0.1);

    SolverParams params;
    params.tolerance = 1e-8;
    params.max_iterations = 500;
    params.verbose = false;
    params.lambda = 0.1;

    int N = 30;
    Array2D rhs(N, N, 1.0);
    Array2D solution(N, N, 0.0);

    // First solve
    solver.solve(rhs, solution, params);
    int iterations1 = solver.get_stats().iterations;

    // Reset
    solver.reset();

    // Second solve (should have fresh statistics)
    solution.fill(0.0);
    solver.solve(rhs, solution, params);
    int iterations2 = solver.get_stats().iterations;

    // Should have same number of iterations
    EXPECT_EQ(iterations1, iterations2);
}

/**
 * @test SmallProblemTest - Test with very small problem
 */
TEST(CGTest, SmallProblemTest) {
    int N = 5;
    double lambda = 0.1;

    Array2D rhs(N, N, 1.0);
    Array2D solution(N, N, 0.0);

    ConjugateGradientSolver solver(false, lambda);
    SolverParams params;
    params.tolerance = 1e-10;
    params.max_iterations = 100;
    params.verbose = false;
    params.lambda = lambda;

    solver.solve(rhs, solution, params);
    const auto& stats = solver.get_stats();

    EXPECT_LT(stats.residual, params.tolerance);
}

/**
 * @test LargeProblemTest - Test with large problem
 */
TEST(CGTest, LargeProblemTest) {
    int N = 150;
    double lambda = 0.05;

    Array2D rhs(N, N, 1.0);
    Array2D solution(N, N, 0.0);

    ConjugateGradientSolver solver(false, lambda);
    SolverParams params;
    params.tolerance = 1e-6;
    params.max_iterations = 2000;
    params.verbose = false;
    params.lambda = lambda;

    auto start = std::chrono::high_resolution_clock::now();
    solver.solve(rhs, solution, params);
    auto end = std::chrono::high_resolution_clock::now();

    const auto& stats = solver.get_stats();
    double time = std::chrono::duration<double>(end - start).count();

    std::cout << "\nLarge problem (N=" << N << "):" << std::endl;
    std::cout << "Iterations: " << stats.iterations << std::endl;
    std::cout << "Time: " << time << " s" << std::endl;
    std::cout << "Time per iteration: " << time / stats.iterations << " s" << std::endl;

    EXPECT_LT(stats.residual, params.tolerance);
}

/**
 * @test LambdaVariationTest - Test with different lambda values
 */
TEST(CGTest, LambdaVariationTest) {
    int N = 40;
    std::vector<double> lambdas = {0.001, 0.01, 0.1, 1.0, 10.0};

    std::cout << "\nLambda variation test:" << std::endl;
    std::cout << "Lambda\tIterations\tResidual\tTime (s)" << std::endl;

    for (double lambda : lambdas) {
        Array2D rhs(N, N, 1.0);
        Array2D solution(N, N, 0.0);

        ConjugateGradientSolver solver(false, lambda);
        SolverParams params;
        params.tolerance = 1e-6;
        params.max_iterations = 1000;
        params.verbose = false;
        params.lambda = lambda;

        auto start = std::chrono::high_resolution_clock::now();
        solver.solve(rhs, solution, params);
        auto end = std::chrono::high_resolution_clock::now();
        double time = std::chrono::duration<double>(end - start).count();

        const auto& stats = solver.get_stats();

        std::cout << std::fixed << lambda << "\t"
                  << stats.iterations << "\t"
                  << std::scientific << stats.residual << "\t"
                  << std::fixed << time << std::endl;

        EXPECT_LT(stats.residual, params.tolerance);
    }
}

// Main function for Google Test
