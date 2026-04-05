/**
 * @file sor_solver.hpp
 * @brief Successive Over-Relaxation (SOR) solver implementation
 *
 * This file implements the SOR iterative method for solving linear systems
 * arising from discretized PDEs. SOR is an accelerated version of Gauss-Seidel
 * with a relaxation parameter omega.
 *
 * The solver supports both serial and parallel (MPI) execution through
 * compile-time template specialization.
 */

#ifndef SOR_SOLVER_HPP
#define SOR_SOLVER_HPP

#include "solver_interface.hpp"
#include <mpi.h>
#include <cmath>
#include <stdexcept>
#include <memory>
#include <chrono>

/**
 * @class SORSolver
 * @brief Successive Over-Relaxation iterative solver
 *
 * SOR is an iterative method for solving linear systems of the form Ax = b.
 * It accelerates the Gauss-Seidel method using a relaxation parameter omega.
 *
 * For omega = 1, SOR reduces to Gauss-Seidel.
 * For 0 < omega < 1, it becomes under-relaxation.
 * For 1 < omega < 2, it becomes over-relaxation (typically faster).
 *
 * The update formula for each grid point is:
 * x[i][j] = (1 - omega) * x_old[i][j] + omega * coeff * (
 *     rhs[i][j] + lambda * (
 *         x[i][j+1] + x[i][j-1] +
 *         x[i+1][j] + x[i-1][j]
 *     )
 * )
 *
 * The solver can use either:
 * - Dictionary (lexicographic) ordering: sequential grid traversal
 * - Red-black ordering: alternates updating red and black points (enables parallelism)
 */
class SORSolver : public ISolver {
public:
    /**
     * @brief Serial constructor
     * @param omega Relaxation parameter (0.0 for automatic optimal calculation)
     * @throws std::invalid_argument if omega >= 2.0 or omega < 0.0
     *
     * omega = 0.0: Automatically compute optimal omega based on grid size
     * 0.0 < omega < 1.0: Under-relaxation (more stable, slower convergence)
     * omega = 1.0: Gauss-Seidel (standard convergence)
     * 1.0 < omega < 2.0: Over-relaxation (faster convergence)
     * omega >= 2.0: Invalid (will diverge)
     */
    explicit SORSolver(double omega = 0.0);

    /**
     * @brief Parallel constructor with MPI support
     * @param omega Relaxation parameter (0.0 for automatic optimal calculation)
     * @param ghost_exchange Pointer to GhostCellExchange for MPI communication
     * @param comm MPI communicator (default: MPI_COMM_WORLD)
     * @throws std::invalid_argument if omega >= 2.0 or omega < 0.0
     *
     * This constructor enables parallel SOR with ghost cell exchange for
     * domain decomposition. The ghost_exchange parameter must be valid and
     * will be used to exchange boundary data between processes.
     */
    SORSolver(double omega, GhostCellExchange* ghost_exchange, MPI_Comm comm = MPI_COMM_WORLD);

    /**
     * @brief Destructor
     */
    ~SORSolver() override = default;

    // Delete copy constructor and copy assignment (solver has MPI state)
    SORSolver(const SORSolver&) = delete;
    SORSolver& operator=(const SORSolver&) = delete;

    // Move constructor and move assignment
    SORSolver(SORSolver&& other) noexcept;
    SORSolver& operator=(SORSolver&& other) noexcept;

    /**
     * @brief Solve the linear system using SOR iteration
     * @param rhs Right-hand side array
     * @param solution Solution array (used as initial guess)
     * @param params Solver parameters
     * @throws std::invalid_argument if array dimensions don't match
     * @throws std::runtime_error if solver encounters an error
     */
    void solve(const utils::Array2D& rhs, utils::Array2D& solution,
               const SolverParams& params) override;

    /**
     * @brief Get solver statistics
     * @return SolverStats from last solve operation
     */
    SolverStats get_stats() const override;

    /**
     * @brief Get solver name
     * @return "SOR" or "SOR-Parallel" or "SOR-RedBlack"
     */
    std::string get_name() const override;

    /**
     * @brief Get solver type
     * @return SolverType::SOR
     */
    SolverType get_type() const override;

    /**
     * @brief Reset solver state
     * Clears cached data and statistics
     */
    void reset() override;

    /**
     * @brief Set relaxation parameter omega
     * @param omega New omega value (0.0 for auto-calculation)
     * @throws std::invalid_argument if omega >= 2.0 or omega < 0.0
     *
     * Setting omega = 0.0 will trigger automatic calculation on next solve.
     */
    void set_omega(double omega);

    /**
     * @brief Get current relaxation parameter
     * @return Current omega value
     */
    double get_omega() const;

    /**
     * @brief Check if omega was automatically computed
     * @return true if omega was computed, false if manually set
     */
    bool is_auto_omega() const;

    /**
     * @brief Enable or disable red-black ordering
     * @param enable true to enable red-black ordering, false for dictionary order
     *
     * Red-black ordering can improve parallelization by allowing
     * independent updates of red and black points.
     */
    void enable_red_black(bool enable);

    /**
     * @brief Check if red-black ordering is enabled
     * @return true if using red-black ordering
     */
    bool is_red_black() const;

private:
    double omega_;                           ///< Relaxation parameter
    double auto_omega_;                     ///< Automatically computed omega
    bool use_auto_omega_;                    ///< Whether to use automatic omega
    bool use_red_black_;                     ///< Whether to use red-black ordering
    bool is_parallel_;                       ///< Whether running in parallel mode
    GhostCellExchange* ghost_exchange_;       ///< Ghost cell exchange for MPI
    MPI_Comm comm_;                          ///< MPI communicator
    int mpi_rank_;                           ///< MPI rank
    int mpi_size_;                           ///< MPI size
    SolverStats stats_;                      ///< Solver statistics

    /**
     * @brief Compute optimal omega for given grid size
     * @param nx Number of grid points in x-direction
     * @param ny Number of grid points in y-direction
     * @return Optimal omega value
     *
     * For the discrete Poisson equation, the optimal relaxation parameter is:
     * omega_opt = 2 / (1 + sqrt(1 - rho^2))
     * where rho = cos(pi/(nx+1)) + cos(pi/(ny+1))
     */
    double compute_optimal_omega(size_t nx, size_t ny);

    /**
     * @brief Validate omega parameter
     * @param omega Omega value to validate
     * @throws std::invalid_argument if omega is invalid
     */
    void validate_omega(double omega);

    /**
     * @brief Compute residual norm ||Ax - b||
     * @param rhs Right-hand side array
     * @param solution Current solution
     * @param lambda Lambda coefficient from discretization
     * @return L-infinity norm of residual
     */
    double compute_residual(const utils::Array2D& rhs,
                            const utils::Array2D& solution,
                            double lambda);

    /**
     * @brief Perform one SOR iteration with dictionary ordering
     * @param rhs Right-hand side array
     * @param solution Solution array (updated in-place)
     * @param lambda Lambda coefficient
     *
     * Updates solution in lexicographic (dictionary) order:
     * (0,0), (0,1), ..., (0,nx-1), (1,0), (1,1), ..., (ny-1,nx-1)
     */
    void sor_iteration_dict(const utils::Array2D& rhs,
                            utils::Array2D& solution,
                            double lambda);

    /**
     * @brief Perform one SOR iteration with red-black ordering
     * @param rhs Right-hand side array
     * @param solution Solution array (updated in-place)
     * @param lambda Lambda coefficient
     * @param color 0 for red points, 1 for black points
     *
     * Updates only points of the specified color:
     * Red:   (i + j) % 2 == 0
     * Black: (i + j) % 2 == 1
     *
     * Red and black points can be updated independently, enabling parallelization.
     */
    void sor_iteration_redblack(const utils::Array2D& rhs,
                                utils::Array2D& solution,
                                double lambda,
                                int color);

    /**
     * @brief Perform parallel SOR iteration with ghost cell exchange
     * @param rhs Right-hand side array
     * @param solution Solution array (updated in-place)
     * @param lambda Lambda coefficient
     * @throws std::runtime_error if ghost_exchange is null
     */
    void sor_iteration_parallel(const utils::Array2D& rhs,
                                utils::Array2D& solution,
                                double lambda);

    /**
     * @brief Validate array dimensions
     * @param rhs Right-hand side array
     * @param solution Solution array
     * @throws std::invalid_argument if dimensions don't match
     */
    void validate_arrays(const utils::Array2D& rhs,
                         const utils::Array2D& solution);

    /**
     * @brief Initialize MPI state
     */
    void init_mpi();

    /**
     * @brief Check if running in parallel
     * @return true if MPI is initialized and size > 1
     */
    bool check_parallel() const;
};

#endif // SOR_SOLVER_HPP
