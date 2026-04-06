/**
 * @file jacobi_solver.hpp
 * @brief Jacobi iterative solver for 2D heat equation
 *
 * This header implements the Jacobi iterative method for solving the linear
 * system arising from the implicit discretization of the 2D heat equation.
 *
 * Features:
 * - Template-based serial/parallel versions
 * - Efficient ghost cell exchange for MPI
 * - Residual computation with global reduction
 * - Performance statistics tracking
 */

#ifndef JACOBI_SOLVER_HPP
#define JACOBI_SOLVER_HPP

#include "solver_interface.hpp"
#include "../../utils/array2d.hpp"
#include "../../mpi/ghost_cell_exchange.hpp"
#include <mpi.h>
#include <memory>
#include <chrono>
#include <cmath>

/**
 * @class JacobiSolver
 * @brief Jacobi solver implementation
 *
 * Implements the Jacobi iterative method for solving (I - lambda*L) * x = b,
 * where L is the discrete Laplacian operator.
 *
 * Jacobi update formula:
 * x_new[i][j] = coeff * (b[i][j] + lambda * (
 *     x_old[i+1][j] + x_old[i][j+1] +
 *     x_old[i-1][j] + x_old[i][j-1]
 * ))
 *
 * where coeff = 1 / (1 + 4*lambda)
 *
 * Template parameter:
 * - Parallel = false: Serial solver (no MPI communication)
 * - Parallel = true:  Parallel solver (with ghost cell exchange)
 *
 * Usage example:
 * @code
 * // Serial version
 * JacobiSolver<false> solver_serial;
 * solver_serial.solve(rhs, solution, params);
 *
 * // Parallel version
 * GhostCellExchange ghost_exchange(nx, ny, topology);
 * JacobiSolver<true> solver_parallel(&ghost_exchange, MPI_COMM_WORLD);
 * solver_parallel.solve(rhs, solution, params);
 * @endcode
 */
template<bool Parallel = false>
class JacobiSolver : public ISolver {
public:
    /**
     * @brief Constructor for serial solver
     *
     * Only available when Parallel = false.
     */
    explicit JacobiSolver() : JacobiSolver(nullptr, MPI_COMM_NULL) {}

    /**
     * @brief Constructor for parallel solver
     * @param ghost_exchange Pointer to ghost cell exchange manager
     * @param comm MPI communicator for parallel operations
     *
     * Only available when Parallel = true.
     */
    JacobiSolver(GhostCellExchange* ghost_exchange, MPI_Comm comm)
        : ghost_exchange_(ghost_exchange)
        , comm_(comm)
        , stats_()
    {
        if constexpr (Parallel) {
            if (ghost_exchange_ == nullptr) {
                throw std::invalid_argument(
                    "Parallel JacobiSolver requires non-null ghost_exchange");
            }
        }
    }

    /**
     * @brief Destructor
     */
    ~JacobiSolver() override = default;

    // Delete copy constructor and copy assignment
    JacobiSolver(const JacobiSolver&) = delete;
    JacobiSolver& operator=(const JacobiSolver&) = delete;

    /**
     * @brief Move constructor
     * @param other Solver to move from
     */
    JacobiSolver(JacobiSolver&& other) noexcept
        : ghost_exchange_(other.ghost_exchange_)
        , comm_(other.comm_)
        , stats_(other.stats_)
    {
        other.ghost_exchange_ = nullptr;
        other.comm_ = MPI_COMM_NULL;
    }

    /**
     * @brief Move assignment operator
     * @param other Solver to move from
     * @return Reference to this solver
     */
    JacobiSolver& operator=(JacobiSolver&& other) noexcept {
        if (this != &other) {
            ghost_exchange_ = other.ghost_exchange_;
            comm_ = other.comm_;
            stats_ = other.stats_;
            other.ghost_exchange_ = nullptr;
            other.comm_ = MPI_COMM_NULL;
        }
        return *this;
    }

    /**
     * @brief Solve the linear system using Jacobi iteration
     * @param rhs Right-hand side array b
     * @param solution Solution array x (will be overwritten)
     * @param params Solver parameters
     *
     * @throws std::invalid_argument if array dimensions don't match
     * @throws ConvergenceFailureException if max iterations exceeded
     */
    void solve(const utils::Array2D& rhs,
               utils::Array2D& solution,
               const SolverParams& params) override {
        auto start_time = std::chrono::high_resolution_clock::now();

        // Validate dimensions
        if (rhs.rows() != solution.rows() || rhs.cols() != solution.cols()) {
            throw std::invalid_argument(
                "rhs and solution must have matching dimensions");
        }

        // Initialize statistics
        stats_.reset();
        stats_.iterations = 0;
        stats_.converged = false;

        // Determine interior dimensions (excluding ghost cells)
        // Assumes arrays include ghost cells: rows = Ny+2, cols = Nx+2
        int Ny = rhs.rows() - 2;  // Interior rows
        int Nx = rhs.cols() - 2;  // Interior columns

        if (Ny <= 0 || Nx <= 0) {
            throw std::invalid_argument(
                "Arrays must be large enough for ghost cells (at least 2x2)");
        }

        // Initialize solution with zeros
        solution.fill(0.0);

        // Working array for previous iteration (x0)
        utils::Array2D x0(rhs.rows(), rhs.cols());
        x0.copy_from(solution);

        // Jacobi coefficient
        double coeff = 1.0 / (1.0 + 4.0 * params.lambda);
        double residual = 1.0;  // Initial residual

        // Main iteration loop
        while (stats_.iterations < params.max_iterations && residual > params.tolerance) {
            stats_.iterations++;

            // Perform one Jacobi iteration
            residual = jacobi_iteration(solution, x0, rhs, coeff, Nx, Ny);

            // Parallel: exchange ghost cells
            if constexpr (Parallel) {
                ghost_exchange_->exchange(solution);
            }

            // Copy new solution to old for next iteration
            x0.copy_from(solution);

            // Parallel: compute global residual
            if constexpr (Parallel) {
                double global_residual;
                MPI_Allreduce(&residual, &global_residual, 1, MPI_DOUBLE, MPI_MAX, comm_);
                residual = global_residual;
            }

            // Check residual periodically if needed
            if (stats_.iterations % params.residual_check_interval == 0) {
                // Could store residual history here if needed
            }
        }

        // Check convergence
        stats_.final_residual = residual;
        stats_.converged = (residual <= params.tolerance);

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        stats_.solve_time = elapsed.count();

        if (!stats_.converged) {
            throw ConvergenceFailureException(
                "Jacobi solver failed to converge after " +
                std::to_string(stats_.iterations) + " iterations. " +
                "Final residual: " + std::to_string(stats_.final_residual) +
                ", tolerance: " + std::to_string(params.tolerance));
        }
    }

    /**
     * @brief Get solver statistics
     * @return SolverStats structure
     */
    SolverStats get_stats() const override {
        return stats_;
    }

    /**
     * @brief Get solver name
     * @return "Jacobi" or "Jacobi-Parallel"
     */
    std::string get_name() const override {
        if constexpr (Parallel) {
            return "Jacobi-Parallel";
        } else {
            return "Jacobi";
        }
    }

    /**
     * @brief Get solver type
     * @return SolverType::Jacobi
     */
    SolverType get_type() const override {
        return SolverType::Jacobi;
    }

    /**
     * @brief Reset solver state
     */
    void reset() override {
        stats_.reset();
    }

private:
    GhostCellExchange* ghost_exchange_;  ///< Ghost cell exchange (parallel only)
    MPI_Comm comm_;                      ///< MPI communicator
    SolverStats stats_;                  ///< Solver statistics

    /**
     * @brief Perform one Jacobi iteration
     * @param x Solution array (output: new values)
     * @param x0 Previous solution array (input)
     * @param b Right-hand side array
     * @param coeff Jacobi coefficient = 1/(1+4*lambda)
     * @param Nx Interior columns
     * @param Ny Interior rows
     * @return Maximum residual for this iteration
     *
     * Updates interior points using the Jacobi formula:
     * x[j][i] = coeff * (lambda * (x0[j+1][i] + x0[j][i+1] +
     *                             x0[j-1][i] + x0[j][i-1]) + b[j][i])
     *
     * Note: Array indices follow row-major convention where first index is row (j)
     * and second index is column (i). This matches the legacy code layout.
     */
    double jacobi_iteration(utils::Array2D& x,
                            const utils::Array2D& x0,
                            const utils::Array2D& b,
                            double coeff,
                            int Nx,
                            int Ny) {
        double max_residual = 0.0;

        // Compute lambda from coeff: coeff = 1/(1+4*lambda)
        double lambda = (1.0 / coeff - 1.0) / 4.0;

        // Loop over interior points (1..Ny for rows, 1..Nx for columns)
        // Note: Ghost cells are at indices 0 and Ny+1 / Nx+1
        for (int j = 1; j <= Ny; ++j) {
            for (int i = 1; i <= Nx; ++i) {
                // Jacobi update formula from legacy code:
                // x[j][i] = Coeff * ( lambda * (x0[j+1][i] + x0[j][i+1] +
                //                               x0[j-1][i] + x0[j][i-1]) + b[j][i] );
                // where Coeff = 1.0/(1.0+4.0*lambda)

                // This solves: (I - lambda*L) * x = b
                // which rearranges to: x_i = (b_i + lambda * sum(neighbors)) / (1 + 4*lambda)

                double new_val = coeff * (
                    lambda * (
                        x0(j+1, i) + x0(j, i+1) +
                        x0(j-1, i) + x0(j, i-1)
                    ) + b(j, i)
                );

                x(j, i) = new_val;

                // Compute residual (max absolute change)
                double residual = std::abs(x(j, i) - x0(j, i));
                max_residual = std::max(max_residual, residual);
            }
        }

        return max_residual;
    }

    /**
     * @brief Compute residual between two arrays
     * @param
     * @param x1 First array
     * @param x2 Second array
     * @return L-infinity norm of difference
     *
     * This method computes max |x1 - x2| over all interior points.
     * Used for convergence checking.
     */
    double compute_residual(const utils::Array2D& x1,
                           const utils::Array2D& x2,
                           int Nx,
                           int Ny) {
        double max_diff = 0.0;

        for (int j = 1; j <= Ny; ++j) {
            for (int i = 1; i <= Nx; ++i) {
                double diff = std::abs(x1(j, i) - x2(j, i));
                max_diff = std::max(max_diff, diff);
            }
        }

        return max_diff;
    }
};

#endif // JACOBI_SOLVER_HPP
