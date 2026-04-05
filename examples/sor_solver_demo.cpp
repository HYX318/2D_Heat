/**
 * @file sor_solver_demo.cpp
 * @brief Demonstration of SOR solver usage
 *
 * This example shows how to use the SOR solver to solve the 2D heat equation.
 */

#include "core/solver/sor_solver.hpp"
#include "utils/array2d.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

/**
 * @brief Create a test problem: heat equation with source term
 */
void create_heat_equation_problem(size_t nx, size_t ny,
                                  utils::Array2D& rhs,
                                  utils::Array2D& solution,
                                  double lambda = 0.25) {
    rhs = utils::Array2D(ny, nx);
    solution = utils::Array2D(ny, nx);
    
    // Initialize solution to zero (initial guess)
    solution.fill(0.0);
    
    // Create a heat source in the center of the domain
    double cx = nx / 2.0;
    double cy = ny / 2.0;
    double sigma = std::min(nx, ny) / 6.0;
    
    for (size_t i = 0; i < ny; ++i) {
        for (size_t j = 0; j < nx; ++j) {
            // Gaussian heat source
            double dx = j - cx;
            double dy = i - cy;
            double r2 = dx*dx + dy*dy;
            rhs(i, j) = 100.0 * std::exp(-r2 / (2.0 * sigma * sigma)) * lambda;
        }
    }
    
    // Set boundary conditions (Dirichlet: u = 0 on boundaries)
    for (size_t j = 0; j < nx; ++j) {
        solution(0, j) = 0.0;
        solution(ny-1, j) = 0.0;
    }
    for (size_t i = 0; i < ny; ++i) {
        solution(i, 0) = 0.0;
        solution(i, nx-1) = 0.0;
    }
}

/**
 * @brief Print solution array
 */
void print_solution(const utils::Array2D& solution, const std::string& title) {
    std::cout << "\n" << title << ":\n";
    std::cout << std::fixed << std::setprecision(3);
    
    for (size_t i = 0; i < solution.rows(); ++i) {
        for (size_t j = 0; j < solution.cols(); ++j) {
            std::cout << std::setw(8) << solution(i, j) << " ";
        }
        std::cout << "\n";
    }
}

/**
 * @brief Print solver statistics
 */
void print_statistics(const SolverStats& stats) {
    std::cout << "\n=== Solver Statistics ===\n";
    std::cout << "Converged:       " << (stats.converged ? "Yes" : "No") << "\n";
    std::cout << "Iterations:      " << stats.iterations << "\n";
    std::cout << "Initial residual: " << std::scientific << stats.initial_residual << "\n";
    std::cout << "Final residual:   " << stats.final_residual << "\n";
    std::cout << "Reduction factor: " << std::fixed << stats.reduction_factor << "\n";
    std::cout << "Solve time:      " << stats.solve_time << " seconds\n";
}

int main(int argc, char** argv) {
    std::cout << "=== SOR Solver Demonstration ===\n\n";
    
    // Problem size
    const size_t nx = 10;
    const size_t ny = 10;
    
    // Solver parameters
    SolverParams params;
    params.tolerance = 1e-6;
    params.max_iterations = 1000;
    params.residual_check_interval = 10;
    params.lambda = 0.25;  // dt / h^2
    
    // Create test problem
    utils::Array2D rhs, solution;
    create_heat_equation_problem(nx, ny, rhs, solution, params.lambda);
    
    std::cout << "Problem size: " << nx << " x " << ny << "\n";
    std::cout << "Tolerance: " << params.tolerance << "\n";
    std::cout << "Max iterations: " << params.max_iterations << "\n";
    
    // Test 1: Gauss-Seidel (omega = 1.0)
    std::cout << "\n--- Test 1: Gauss-Seidel (omega = 1.0) ---\n";
    {
        utils::Array2D sol_copy(solution);
        SORSolver solver(1.0);
        
        std::cout << "Solver name: " << solver.get_name() << "\n";
        std::cout << "Omega: " << solver.get_omega() << "\n";
        
        solver.solve(rhs, sol_copy, params);
        print_statistics(solver.get_stats());
    }
    
    // Test 2: SOR with omega = 1.5
    std::cout << "\n--- Test 2: SOR with omega = 1.5 ---\n";
    {
        utils::Array2D sol_copy(solution);
        SORSolver solver(1.5);
        
        std::cout << "Solver name: " << solver.get_name() << "\n";
        std::cout << "Omega: " << solver.get_omega() << "\n";
        
        solver.solve(rhs, sol_copy, params);
        print_statistics(solver.get_stats());
    }
    
    // Test 3: SOR with red-black ordering
    std::cout << "\n--- Test 3: SOR with red-black ordering (omega = 1.5) ---\n";
    {
        utils::Array2D sol_copy(solution);
        SORSolver solver(1.5);
        solver.enable_red_black(true);
        
        std::cout << "Solver name: " << solver.get_name() << "\n";
        std::cout << "Omega: " << solver.get_omega() << "\n";
        std::cout << "Red-black ordering: " << (solver.is_red_black() ? "enabled" : "disabled") << "\n";
        
        solver.solve(rhs, sol_copy, params);
        print_statistics(solver.get_stats());
    }
    
    // Test 4: SOR with automatic omega
    std::cout << "\n--- Test 4: SOR with automatic omega ---\n";
    {
        utils::Array2D sol_copy(solution);
        SORSolver solver(0.0);  // Auto omega
        
        std::cout << "Solver name: " << solver.get_name() << "\n";
        std::cout << "Auto omega: " << (solver.is_auto_omega() ? "yes" : "no") << "\n";
        
        solver.solve(rhs, sol_copy, params);
        
        std::cout << "Computed omega: " << solver.get_omega() << "\n";
        print_statistics(solver.get_stats());
    }
    
    // Test 5: Compare different omega values
    std::cout << "\n--- Test 5: Comparison of different omega values ---\n";
    std::cout << std::setw(10) << "Omega" 
              << std::setw(12) << "Iterations" 
              << std::setw(15) << "Final Residual"
              << std::setw(12) << "Time (s)\n";
    std::cout << std::string(50, '-') << "\n";
    
    double omegas[] = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8};
    for (double omega : omegas) {
        utils::Array2D sol_copy(solution);
        SORSolver solver(omega);
        solver.solve(rhs, sol_copy, params);
        
        SolverStats stats = solver.get_stats();
        std::cout << std::fixed << std::setprecision(2);
        std::cout << std::setw(10) << omega
                  << std::setw(12) << stats.iterations
                  << std::scientific << std::setprecision(2)
                  << std::setw(15) << stats.final_residual
                  << std::fixed << std::setprecision(4)
                  << std::setw(12) << stats.solve_time << "\n";
    }
    
    std::cout << "\n=== Demonstration Complete ===\n";
    
    return 0;
}
