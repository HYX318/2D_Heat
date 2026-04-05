/**
 * @file test_jacobi_solver.cpp
 * @brief Unit tests for Jacobi solver
 *
 * Test suite for Jacobi solver implementation.
 */

#include <gtest/gtest.h>
#include <mpi.h>
#include <cmath>
#include <memory>

#include "src/core/solver/solver_interface.hpp"
#include "src/core/solver/jacobi_solver.hpp"
#include "src/mpi/ghost_cell_exchange.hpp"
#include "src/mpi/cartesian_topology.hpp"
#include "src/utils/array2d.hpp"

using namespace utils;

/**
 * Test fixture for Jacobi solver tests
 */
class JacobiSolverTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize MPI if not already initialized
        int flag;
        MPI_Initialized(&flag);
        if (!flag) {
            MPI_Init(nullptr, nullptr);
            mpi_initialized_ = true;
        }

        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        MPI_Comm_size(MPI_COMM_WORLD, &size_);
    }

    void TearDown() override {
        // Only finalize MPI if we initialized it
        if (mpi_initialized_) {
            MPI_Finalize();
        }
    }

    /**
     * Create a simple test problem with known solution
     * @param Nx Number of interior columns
     * @param Ny Number of interior rows
     * @param lambda Coefficient parameter
     * @return Pair of (rhs, exact_solution)
     */
    std::pair<Array2D, Array2D> create_test_problem(int Nx, int Ny, double lambda) {
        // Arrays with ghost cells: Ny+2 rows, Nx+2 cols
        Array2D rhs(Ny + 2, Nx + 2);
        Array2D exact(Ny + 2, Nx + 2);

        // Create a simple problem with sinusoidal solution
        double pi = 3.14159265358979323846;

        for (int j = 0; j < Ny + 2; ++j) {
            for (int i = 0; i < Nx + 2; ++i) {
                double x = static_cast<double>(i) / (Nx + 1);
                double y = static_cast<double>(j) / (Ny + 1);

                exact(j, i) = std::sin(pi * x) * std::sin(pi * y);

                // For interior points, compute RHS
                // (I - λ*L)u = b
                if (i >= 1 && i <= Nx && j >= 1 && j <= Ny) {
                    double laplacian = exact(j+1, i) - 2.0*exact(j, i) + exact(j-1, i) +
                                    exact(j, i+1) - 2.0*exact(j, i) + exact(j, i-1);
                    rhs(j, i) = exact(j, i) - lambda * laplacian;
                } else {
                    // Boundary
                    rhs(j, i) = exact(j, i);
                }
            }
        }

        return {rhs, exact};
    }

    /**
     * Compute L-infinity norm of error
     */
    double compute_max_error(const Array2D& computed, const Array2D& exact) {
        double max_error = 0.0;
        for (size_t j = 0; j < computed.rows(); ++j) {
            for (size_t i = 0; i < computed.cols(); ++i) {
                double error = std::abs(computed(j, i) - exact(j, i));
                max_error = std::max(max_error, error);
            }
        }
        return max_error;
    }

    int rank_ = 0;
    int size_ = 1;
    bool mpi_initialized_ = false;
};

/**
 * Test: Convergence behavior
 */
TEST_F(JacobiSolverTest, ConvergenceTest) {
    // Only run on rank 0 for serial tests
    if (rank_ != 0) return;

    int Nx = 10;
    int Ny = 10;
    double lambda = 0.1;

    auto [rhs, exact] = create_test_problem(Nx, Ny, lambda);

    // Create solver
    JacobiSolver<false> solver;

    // Create solution array
    Array2D solution(Ny + 2, Nx + 2);

    // Set parameters
    SolverParams params;
    params.tolerance = 1e-8;
    params.max_iterations = 10000;
    params.lambda = lambda;

    // Solve
    solver.solve(rhs, solution, params);

    // Check convergence
    SolverStats stats = solver.get_stats();
    EXPECT_TRUE(stats.converged);
    EXPECT_LT(stats.final_residual, params.tolerance);

    // Check accuracy
    double max_error = compute_max_error(solution, exact);
    EXPECT_LT(max_error, 1e-6) << "Max error: " << max_error;
}

/**
 * Test: Serial solver basic functionality
 */
TEST_F(JacobiSolverTest, SerialSolverBasicTest) {
    if (rank_ != 0) return;

    int Nx = 5;
    int Ny = 5;
    double lambda = 0.25;

    auto [rhs, exact] = create_test_problem(Nx, Ny, lambda);

    JacobiSolver<false> solver;
    Array2D solution(Ny + 2, Nx + 2);

    SolverParams params;
    params.tolerance = 1e-6;
    params.max_iterations = 1000;
    params.lambda = lambda;

    solver.solve(rhs, solution, params);

    SolverStats stats = solver.get_stats();
    EXPECT_TRUE(stats.converged);
    EXPECT_GT(stats.iterations, 0);
    EXPECT_GT(stats.solve_time, 0.0);

    // Check solver name and type
    EXPECT_EQ(solver.get_name(), "Jacobi");
    EXPECT_EQ(solver.get_type(), SolverType::Jacobi);
}

/**
 * Test: Reset functionality
 */
TEST_F(JacobiSolverTest, ResetTest) {
    if (rank_ != 0) return;

    int Nx = 6;
    int Ny = 6;
    double lambda = 0.15;

    auto [rhs, exact] = create_test_problem(Nx, Ny, lambda);

    JacobiSolver<false> solver;
    Array2D solution(Ny + 2, Nx + 2);

    SolverParams params;
    params.tolerance = 1e-6;
    params.max_iterations = 1000;
    params.lambda = lambda;

    // First solve
    solver.solve(rhs, solution, params);
    auto stats1 = solver.get_stats();

    // Reset and solve again
    solver.reset();
    solution.fill(0.0);
    solver.solve(rhs, solution, params);
    auto stats2 = solver.get_stats();

    // Should produce same results
    EXPECT_TRUE(stats2.converged);
    EXPECT_EQ(stats1.iterations, stats2.iterations);
}

/**
 * Test: Parallel solver initialization
 */
TEST_F(JacobiSolverTest, ParallelSolverInitializationTest) {
    // Create Cartesian topology
    CartesianTopology topology(MPI_COMM_WORLD);

    // Create domain decomposition (small domain for testing)
    int nx = 4;  // Interior columns per process
    int ny = 4;  // Interior rows per process

    GhostCellExchange ghost_exchange(nx, ny, topology);

    // Create solver (should succeed on all ranks)
    JacobiSolver<true> solver(&ghost_exchange, MPI_COMM_WORLD);

    // Check solver properties
    EXPECT_EQ(solver.get_name(), "Jacobi-Parallel");
    EXPECT_EQ(solver.get_type(), SolverType::Jacobi);
}

/**
 * Test: Residual computation accuracy
 */
TEST_F(JacobiSolverTest, ResidualComputationTest) {
    if (rank_ != 0) return;

    int Nx = 6;
    int Ny = 6;
    double lambda = 0.3;

    auto [rhs, exact] = create_test_problem(Nx, Ny, lambda);

    JacobiSolver<false> solver;
    Array2D solution(Ny + 2, Nx + 2);

    // Solve with a known tolerance
    SolverParams params;
    params.tolerance = 1e-10;
    params.max_iterations = 2000;
    params.lambda = lambda;

    solver.solve(rhs, solution, params);

    SolverStats stats = solver.get_stats();

    // Final residual should be below tolerance
    EXPECT_LT(stats.final_residual, params.tolerance);

    // Compute actual residual by applying operator
    double max_resid = 0.0;
    double coeff = 1.0 / (1.0 + 4.0 * lambda);

    for (int j = 1; j <= Ny; ++j) {
        for (int i = 1; i <= Nx; ++i) {
            double lhs = solution(j, i) - lambda * (
                solution(j+1, i) - 2.0 * solution(j, i) + solution(j-1, i) +
                solution(j, i+1) - 2.0 * solution(j, i) + solution(j, i-1)
            );
            double resid = std::abs(lhs - rhs(j, i));
            max_resid = std::max(max_resid, resid);
        }
    }

    // Computed residual should match solver's reported residual
    EXPECT_NEAR(max_resid, stats.final_residual, 1e-8);
}

/**
 * Test: Edge cases
 */
TEST_F(JacobiSolverTest, EdgeCasesTest) {
    if (rank_ != 0) return;

    // Very small problem
    {
        int Nx = 2;
        int Ny = 2;
        double lambda = 0.1;

        auto [rhs, exact] = create_test_problem(Nx, Ny, lambda);

        JacobiSolver<false> solver;
        Array2D solution(Ny + 2, Nx + 2);

        SolverParams params;
        params.tolerance = 1e-6;
        params.max_iterations = 1000;
        params.lambda = lambda;

        solver.solve(rhs, solution, params);
        EXPECT_TRUE(solver.get_stats().converged);
    }

    // Very small lambda
    {
        int Nx = 5;
        int Ny = 5;
        double lambda = 0.001;

        auto [rhs, exact] = create_test_problem(Nx, Ny, lambda);

        JacobiSolver<false> solver;
        Array2D solution(Ny + 2, Nx + 2);

        SolverParams params;
        params.tolerance = 1e-6;
        params.max_iterations = 1000;
        params.lambda = lambda;

        solver.solve(rhs, solution, params);
        EXPECT_TRUE(solver.get_stats().converged);
    }

    // Larger lambda (but still stable)
    {
        int Nx = 5;
        int Ny = 5;
        double lambda = 0.4;  // Still < 0.5 for stability

        auto [rhs, exact] = create_test_problem(Nx, Ny, lambda);

        JacobiSolver<false> solver;
        Array2D solution(Ny + 2, Nx + 2);

        SolverParams params;
        params.tolerance = 1e-6;
        params.max_iterations = 2000;  // May need more iterations
        params.lambda = lambda;

        solver.solve(rhs, solution, params);
        EXPECT_TRUE(solver.get_stats().converged);
    }
}

/**
 * Test: Performance benchmark
 */
TEST_F(JacobiSolverTest, PerformanceTest) {
    if (rank_ != 0) return;

    int Nx = 50;
    int Ny = 50;
    double lambda = 0.2;

    auto [rhs, exact] = create_test_problem(Nx, Ny, lambda);

    JacobiSolver<false> solver;
    Array2D solution(Ny + 2, Nx + 2);

    SolverParams params;
    params.tolerance = 1e-6;
    params.max_iterations = 10000;
    params.lambda = lambda;

    auto start = std::chrono::high_resolution_clock::now();
    solver.solve(rhs, solution, params);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    double solve_time = elapsed.count();

    SolverStats stats = solver.get_stats();

    std::cout << "Performance Test Results:" << std::endl;
    std::cout << "  Problem size: " << Nx << "x" << Ny << std::endl;
    std::cout << "  Iterations: " << stats.iterations << std::endl;
    std::cout << "  Solve time: " << solve_time << " seconds" << std::endl;
    std::cout << "  Time per iteration: " << solve_time / stats.iterations << " seconds" << std::endl;
    std::cout << "  Final residual: " << stats.final_residual << std::endl;

    // Basic sanity checks
    EXPECT_TRUE(stats.converged);
    EXPECT_GT(stats.iterations, 0);
    EXPECT_GT(solve_time, 0.0);
}

/**
 * Main function for MPI test execution
 */
int main(int argc, char** argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Initialize Google Test
    ::testing::InitGoogleTest(&argc, argv);

    // Run tests
    int result = RUN_ALL_TESTS();

    // Finalize MPI
    MPI_Finalize();

    return result;
}
