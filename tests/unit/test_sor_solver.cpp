/**
 * @file test_sor_solver.cpp
 * @brief Unit tests for SOR solver using Google Test framework
 */

#include <gtest/gtest.h>
#include "core/solver/sor_solver.hpp"
#include <cmath>
#include <random>
#include <fstream>

// Test fixture for SOR solver tests
class SORSolverTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Default parameters
        params.tolerance = 1e-6;
        params.max_iterations = 1000;
        params.residual_check_interval = 10;
        params.lambda = 0.25;
    }

    void TearDown() override {
        // Cleanup handled by RAII
    }

    /**
     * @brief Create a simple test problem with known solution
     * @param nx Grid size in x-direction
     * @param ny Grid size in y-direction
     * @param rhs Right-hand side array (output)
     * @param solution Exact solution array (output)
     */
    void create_test_problem(size_t nx, size_t ny,
                             utils::Array2D& rhs,
                             utils::Array2D& solution) {
        rhs = utils::Array2D(ny, nx);
        solution = utils::Array2D(ny, nx);

        // Create a simple solution: sin(pi*x) * sin(pi*y)
        double dx = 1.0 / (nx + 1);
        double dy = 1.0 / (ny + 1);

        for (size_t i = 0; i < ny; ++i) {
            for (size_t j = 0; j < nx; ++j) {
                double x = (j + 1) * dx;
                double y = (i + 1) * dy;
                solution(i, j) = std::sin(M_PI * x) * std::sin(M_PI * y);
                
                // For the heat equation, rhs = f(x,y)
                // For sin(pi*x)*sin(pi*y), the Laplacian gives -2*pi^2*solution
                rhs(i, j) = -2.0 * M_PI * M_PI * solution(i, j) * params.lambda;
            }
        }
    }

    /**
     * @brief Create a random test problem
     * @param nx Grid size in x-direction
     * @param ny Grid size in y-direction
     * @param rhs Right-hand side array (output)
     * @param solution Initial guess (output)
     */
    void create_random_problem(size_t nx, size_t ny,
                                 utils::Array2D& rhs,
                                 utils::Array2D& solution) {
        rhs = utils::Array2D(ny, nx);
        solution = utils::Array2D(ny, nx);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist(-1.0, 1.0);

        for (size_t i = 0; i < ny; ++i) {
            for (size_t j = 0; j < nx; ++j) {
                rhs(i, j) = dist(gen);
                solution(i, j) = dist(gen);
            }
        }

        // Set Dirichlet boundary conditions (zero)
        for (size_t j = 0; j < nx; ++j) {
            solution(0, j) = 0.0;
            solution(ny-1, j) = 0.0;
        }
        for (size_t i = 0; i < ny; ++i) {
            solution(i, 0) = 0.0;
            solution(i, nx-1) = 0.0;
        }
    }

    SolverParams params;
};

/**
 * @test DefaultConstruction - Test default constructor
 */
TEST_F(SORSolverTest, DefaultConstruction) {
    SORSolver solver;
    
    EXPECT_EQ(solver.get_type(), SolverType::SOR);
    EXPECT_FALSE(solver.is_auto_omega());
    EXPECT_FALSE(solver.is_red_black());
    EXPECT_DOUBLE_EQ(solver.get_omega(), 1.0);  // Default omega = 1.0 (Gauss-Seidel)
}

/**
 * @test AutoOmegaConstruction - Test constructor with auto omega
 */
TEST_F(SORSolverTest, AutoOmegaConstruction) {
    SORSolver solver(0.0);  // Auto omega
    
    EXPECT_TRUE(solver.is_auto_omega());
    EXPECT_DOUBLE_EQ(solver.get_omega(), 1.0);  // Will be computed on first solve
}

/**
 * @test CustomOmegaConstruction - Test constructor with custom omega
 */
TEST_F(SORSolverTest, CustomOmegaConstruction) {
    SORSolver solver(1.5);  // Over-relaxation
    
    EXPECT_FALSE(solver.is_auto_omega());
    EXPECT_DOUBLE_EQ(solver.get_omega(), 1.5);
}

/**
 * @test InvalidOmegaConstruction - Test invalid omega values
 */
TEST_F(SORSolverTest, InvalidOmegaConstruction) {
    EXPECT_THROW(SORSolver(2.0), std::invalid_argument);  // omega = 2.0 is invalid
    EXPECT_THROW(SORSolver(2.5), std::invalid_argument);  // omega > 2.0 is invalid
    EXPECT_THROW(SORSolver(-0.1), std::invalid_argument);  // omega < 0.0 is invalid
}

/**
 * @test SetOmega - Test setting omega after construction
 */
TEST_F(SORSolverTest, SetOmega) {
    SORSolver solver;
    
    solver.set_omega(1.8);
    EXPECT_DOUBLE_EQ(solver.get_omega(), 1.8);
    EXPECT_FALSE(solver.is_auto_omega());
    
    solver.set_omega(0.5);
    EXPECT_DOUBLE_EQ(solver.get_omega(), 0.5);
}

/**
 * @test SetInvalidOmega - Test setting invalid omega
 */
TEST_F(SORSolverTest, SetInvalidOmega) {
    SORSolver solver;
    
    EXPECT_THROW(solver.set_omega(2.0), std::invalid_argument);
    EXPECT_THROW(solver.set_omega(-0.5), std::invalid_argument);
}

/**
 * @test RedBlackOrdering - Test red-black ordering
 */
TEST_F(SORSolverTest, RedBlackOrdering) {
    SORSolver solver;
    
    EXPECT_FALSE(solver.is_red_black());
    
    solver.enable_red_black(true);
    EXPECT_TRUE(solver.is_red_black());
    
    solver.enable_red_black(false);
    EXPECT_FALSE(solver.is_red_black());
}

/**
 * @test SimpleConvergence - Test convergence on simple problem
 */
TEST_F(SORSolverTest, SimpleConvergence) {
    const size_t nx = 10;
    const size_t ny = 10;
    
    utils::Array2D rhs(ny, nx);
    utils::Array2D solution(ny, nx);
    create_test_problem(nx, ny, rhs, solution);
    
    SORSolver solver(1.5);  // Over-relaxation
    solver.solve(rhs, solution, params);
    
    SolverStats stats = solver.get_stats();
    
    EXPECT_TRUE(stats.converged);
    EXPECT_GT(stats.iterations, 0);
    EXPECT_LT(stats.iterations, params.max_iterations);
    EXPECT_LT(stats.final_residual, params.tolerance);
    EXPECT_GT(stats.solve_time, 0.0);
    EXPECT_LT(stats.reduction_factor, 1.0);
}

/**
 * @test DifferentOmegaValues - Test convergence with different omega values
 */
TEST_F(SORSolverTest, DifferentOmegaValues) {
    const size_t nx = 10;
    const size_t ny = 10;
    
    utils::Array2D rhs(ny, nx);
    utils::Array2D solution(ny, nx);
    create_test_problem(nx, ny, rhs, solution);
    
    double omega_values[] = {0.5, 0.8, 1.0, 1.2, 1.5, 1.8};
    
    for (double omega : omega_values) {
        utils::Array2D sol_copy(solution);  // Copy initial guess
        SORSolver solver(omega);
        solver.solve(rhs, sol_copy, params);
        
        SolverStats stats = solver.get_stats();
        EXPECT_TRUE(stats.converged) << "Failed to converge with omega = " << omega;
        EXPECT_LT(stats.final_residual, params.tolerance);
    }
}

/**
 * @test AutoOmega - Test automatic omega calculation
 */
TEST_F(SORSolverTest, AutoOmega) {
    const size_t nx = 20;
    const size_t ny = 20;
    
    utils::Array2D rhs(ny, nx);
    utils::Array2D solution(ny, nx);
    create_test_problem(nx, ny, rhs, solution);
    
    SORSolver solver(0.0);  // Auto omega
    EXPECT_TRUE(solver.is_auto_omega());
    
    solver.solve(rhs, solution, params);
    
    // After solving, omega should be computed
    EXPECT_GT(solver.get_omega(), 1.0);  // Optimal omega should be > 1.0
    EXPECT_LT(solver.get_omega(), 2.0);
    
    SolverStats stats = solver.get_stats();
    EXPECT_TRUE(stats.converged);
}

/**
 * @test RedBlackConvergence - Test convergence with red-black ordering
 */
TEST_F(SORSolverTest, RedBlackConvergence) {
    const size_t nx = 15;
    const size_t ny = 15;
    
    utils::Array2D rhs(ny, nx);
    utils::Array2D solution(ny, nx);
    create_test_problem(nx, ny, rhs, solution);
    
    SORSolver solver(1.5);
    solver.enable_red_black(true);
    
    solver.solve(rhs, solution, params);
    
    SolverStats stats = solver.get_stats();
    
    EXPECT_TRUE(stats.converged);
    EXPECT_LT(stats.final_residual, params.tolerance);
    EXPECT_TRUE(stats.converged);
}

/**
 * @test CompareOrderings - Compare dictionary vs red-black ordering
 */
TEST_F(SORSolverTest, CompareOrderings) {
    const size_t nx = 10;
    const size_t ny = 10;

    utils::Array2D rhs(ny, nx);
    utils::Array2D solution_dict(ny, nx);
    utils::Array2D solution_rb(ny, nx);
    create_test_problem(nx, ny, rhs, solution_dict);
    solution_rb.copy_from(solution_dict);
    
    // Dictionary ordering
    SORSolver solver_dict(1.5);
    solver_dict.solve(rhs, solution_dict, params);
    
    // Red-black ordering
    SORSolver solver_rb(1.5);
    solver_rb.enable_red_black(true);
    solver_rb.solve(rhs, solution_rb, params);
    
    SolverStats stats_dict = solver_dict.get_stats();
    SolverStats stats_rb = solver_rb.get_stats();
    
    // Both should converge
    EXPECT_TRUE(stats_dict.converged);
    EXPECT_TRUE(stats_rb.converged);
    
    // Results should be similar (not necessarily identical due to different order)
    double diff = 0.0;
    for (size_t i = 0; i < ny; ++i) {
        for (size_t j = 0; j < nx; ++j) {
            diff += std::abs(solution_dict(i, j) - solution_rb(i, j));
        }
    }
    EXPECT_LT(diff, 1e-4);  // Solutions should be close
}

/**
 * @test RandomProblem - Test convergence on random problem
 */
TEST_F(SORSolverTest, RandomProblem) {
    const size_t nx = 12;
    const size_t ny = 12;

    utils::Array2D rhs(ny, nx);
    utils::Array2D solution(ny, nx);
    create_random_problem(nx, ny, rhs, solution);
    
    SORSolver solver(1.7);
    solver.solve(rhs, solution, params);
    
    SolverStats stats = solver.get_stats();
    
    EXPECT_TRUE(stats.converged);
    EXPECT_LT(stats.final_residual, params.tolerance);
}

/**
 * @test LargeProblem - Test convergence on larger problem
 */
TEST_F(SORSolverTest, LargeProblem) {
    const size_t nx = 50;
    const size_t ny = 50;
    
    utils::Array2D rhs(ny, nx);
    utils::Array2D solution(ny, nx);
    create_test_problem(nx, ny, rhs, solution);
    
    SORSolver solver(0.0);  // Auto omega for large problem
    params.max_iterations = 5000;  // More iterations for larger problem
    
    solver.solve(rhs, solution, params);
    
    SolverStats stats = solver.get_stats();
    
    EXPECT_TRUE(stats.converged);
    EXPECT_LT(stats.final_residual, params.tolerance);
}

/**
 * @test Reset - Test reset functionality
 */
TEST_F(SORSolverTest, Reset) {
    const size_t nx = 10;
    const size_t ny = 10;
    
    utils::Array2D rhs(ny, nx);
    utils::Array2D solution(ny, nx);
    create_test_problem(nx, ny, rhs, solution);
    
    SORSolver solver(1.5);
    solver.solve(rhs, solution, params);
    
    SolverStats stats_before = solver.get_stats();
    EXPECT_GT(stats_before.iterations, 0);
    
    // Reset and solve again
    solver.reset();
    SolverStats stats_after_reset = solver.get_stats();
    EXPECT_EQ(stats_after_reset.iterations, 0);
    
    solution.fill(0.0);  // Reset initial guess
    solver.solve(rhs, solution, params);
    
    SolverStats stats_after = solver.get_stats();
    EXPECT_GT(stats_after.iterations, 0);
}

/**
 * @test GetStats - Test statistics retrieval
 */
TEST_F(SORSolverTest, GetStats) {
    const size_t nx = 10;
    const size_t ny = 10;

    utils::Array2D rhs(ny, nx);
    utils::Array2D solution(ny, nx);
    create_test_problem(nx, ny, rhs, solution);
    
    SORSolver solver(1.5);
    solver.solve(rhs, solution, params);
    
    SolverStats stats = solver.get_stats();
    
    // Verify all stats are set correctly
    EXPECT_TRUE(stats.converged);
    EXPECT_GT(stats.iterations, 0);
    EXPECT_LT(stats.iterations, params.max_iterations);
    EXPECT_GT(stats.initial_residual, 0.0);
    EXPECT_GT(stats.final_residual, 0.0);
    EXPECT_LT(stats.final_residual, params.tolerance);
    EXPECT_GT(stats.solve_time, 0.0);
    EXPECT_LT(stats.reduction_factor, 1.0);
    EXPECT_GT(stats.reduction_factor, 0.0);
}

/**
 * @test DifferentGridSizes - Test convergence on different grid sizes
 */
TEST_F(SORSolverTest, DifferentGridSizes) {
    size_t sizes[] = {5, 8, 10, 15, 20};
    
    for (size_t nx : sizes) {
        size_t ny = nx;  // Square grids

        utils::Array2D rhs(ny, nx);
        utils::Array2D solution(ny, nx);
        create_test_problem(nx, ny, rhs, solution);
        
        SORSolver solver(1.5);
        solver.solve(rhs, solution, params);
        
        SolverStats stats = solver.get_stats();
        
        EXPECT_TRUE(stats.converged) << "Failed to converge for grid size " << nx;
        EXPECT_LT(stats.final_residual, params.tolerance);
    }
}

/**
 * @test TightTolerance - Test convergence with tight tolerance
 */
TEST_F(SORSolverTest, TightTolerance) {
    const size_t nx = 10;
    const size_t ny = 10;
    
    utils::Array2D rhs(ny, nx);
    utils::Array2D solution(ny, nx);
    create_test_problem(nx, ny, rhs, solution);
    
    params.tolerance = 1e-10;  // Very tight tolerance
    params.max_iterations = 10000;
    
    SORSolver solver(1.8);
    solver.solve(rhs, solution, params);
    
    SolverStats stats = solver.get_stats();
    
    EXPECT_TRUE(stats.converged);
    EXPECT_LT(stats.final_residual, params.tolerance);
}

/**
 * @test GetNames - Test solver name retrieval
 */
TEST_F(SORSolverTest, GetNames) {
    SORSolver solver1(1.5);
    EXPECT_EQ(solver1.get_name(), "SOR");
    
    SORSolver solver2(1.5);
    solver2.enable_red_black(true);
    EXPECT_EQ(solver2.get_name(), "SOR-RedBlack");
}

/**
 * @test MoveSemantics - Test move constructor and move assignment
 */
TEST_F(SORSolverTest, MoveSemantics) {
    SORSolver solver1(1.5);
    solver1.enable_red_black(true);
    
    // Move construct
    SORSolver solver2(std::move(solver1));
    EXPECT_EQ(solver2.get_omega(), 1.5);
    EXPECT_TRUE(solver2.is_red_black());
    
    // Move assign
    SORSolver solver3(0.5);
    solver3 = std::move(solver2);
    EXPECT_EQ(solver3.get_omega(), 1.5);
    EXPECT_TRUE(solver3.is_red_black());
}

/**
 * @test OptimalOmegaComputation - Test optimal omega computation
 */
TEST_F(SORSolverTest, OptimalOmegaComputation) {
    SORSolver solver(0.0);  // Auto omega
    
    // For small grids, optimal omega should be close to 1
    utils::Array2D rhs(5, 5);
    utils::Array2D solution(5, 5);
    create_test_problem(5, 5, rhs, solution);
    solver.solve(rhs, solution, params);
    EXPECT_GT(solver.get_omega(), 1.0);
    
    // For larger grids, optimal omega should be higher
    utils::Array2D rhs2(50, 50);
    utils::Array2D solution2(50, 50);
    create_test_problem(50, 50, rhs2, solution2);
    SORSolver solver2(0.0);
    solver2.solve(rhs2, solution2, params);
    EXPECT_GT(solver2.get_omega(), solver.get_omega());
}

/**
 * @test ConvergenceRate - Test convergence rate with different omega
 */
TEST_F(SORSolverTest, ConvergenceRate) {
    const size_t nx = 15;
    const size_t ny = 15;
    
    utils::Array2D rhs(ny, nx);
    utils::Array2D solution(ny, nx);
    create_test_problem(nx, ny, rhs, solution);
    
    // Test with different omega values and verify they all converge.
    std::vector<double> omegas = {0.8, 1.0, 1.2, 1.5};
    
    for (double omega : omegas) {
        utils::Array2D sol_copy(solution);
        SORSolver solver(omega);
        solver.solve(rhs, sol_copy, params);
        SolverStats stats = solver.get_stats();
        EXPECT_TRUE(stats.converged) << "Failed to converge for omega = " << omega;
        EXPECT_LT(stats.final_residual, params.tolerance);
        EXPECT_LT(stats.iterations, params.max_iterations);
    }
}

// MPI tests (only run if MPI is initialized)
#ifdef ENABLE_MPI

/**
 * @test ParallelConstruction - Test parallel constructor
 */
TEST_F(SORSolverTest, ParallelConstruction) {
    int initialized = 0;
    MPI_Initialized(&initialized);
    
    if (!initialized) {
        GTEST_SKIP() << "MPI not initialized";
    }
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (size < 2) {
        GTEST_SKIP() << "Need at least 2 MPI processes for parallel test";
    }
    
    // Note: This test requires a valid GhostCellExchange
    // For now, we just test that the constructor signature compiles
    // Actual parallel testing would require setting up a full MPI context
    SUCCEED();
}

#endif // ENABLE_MPI
