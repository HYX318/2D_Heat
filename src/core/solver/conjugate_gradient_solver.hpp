/**
 * @file conjugate_gradient_solver.hpp
 * @brief Conjugate Gradient solver for 2D heat equation
 *
 * Implements the Conjugate Gradient (CG) method for solving
 * linear systems arising from implicit Euler discretization
 * of the 2D heat equation.
 *
 * The solver handles:
 * - Standard CG algorithm
 * - Preconditioned CG (PCG) with Jacobi preconditioner
 * - Parallel execution with MPI
 * - Restart mechanism for numerical stability
 */

#ifndef CONJUGATE_GRADIENT_SOLVER_HPP
#define CONJUGATE_GRADIENT_SOLVER_HPP

#include <mpi.h>
#include <memory>
#include <string>
#include "solver_interface.hpp"
#include "../../utils/array2d.hpp"

/**
 * @class ConjugateGradientSolver
 * @brief Conjugate Gradient solver for 2D heat equation
 *
 * Solves the linear system Ax = b where:
 * - A = (I - λ∇²) for implicit Euler discretization
 * - λ = dt / h² is the diffusion coefficient
 *
 * The solver implements the CG algorithm with:
 * - 5-point stencil for matrix-vector multiplication
 * - MPI parallelization for global dot products
 * - Optional Jacobi preconditioning
 * - Restart mechanism for numerical stability
 *
 * Theoretical convergence:
 * - Iterations ~ O(√κ) where κ is condition number
 * - For Poisson: κ = O(N²), so CG needs O(N) iterations
 * - Significantly faster than stationary methods (Jacobi, SOR)
 */
class ConjugateGradientSolver {
public:
    /**
     * @brief Constructor for serial execution
     * @param use_preconditioner Enable Jacobi preconditioner (default: false)
     * @param lambda Diffusion coefficient dt/h² (default: 0.0)
     */
    explicit ConjugateGradientSolver(bool use_preconditioner = false, double lambda = 0.0);

    /**
     * @brief Constructor for parallel execution
     * @param use_preconditioner Enable Jacobi preconditioner
     * @param lambda Diffusion coefficient dt/h²
     * @param comm MPI communicator (default: MPI_COMM_WORLD)
     * @param neighbor_rank Array of neighbor process ranks [N, S, E, W]
     * @param nx Local grid size in x-direction
     * @param ny Local grid size in y-direction
     */
    ConjugateGradientSolver(bool use_preconditioner, double lambda,
                            MPI_Comm comm,
                            const int neighbor_rank[4],
                            int nx, int ny);

    /**
     * @brief Destructor
     */
    ~ConjugateGradientSolver();

    // Delete copy operations
    ConjugateGradientSolver(const ConjugateGradientSolver&) = delete;
    ConjugateGradientSolver& operator=(const ConjugateGradientSolver&) = delete;

    // Allow move operations
    ConjugateGradientSolver(ConjugateGradientSolver&&) noexcept = default;
    ConjugateGradientSolver& operator=(ConjugateGradientSolver&&) noexcept = default;

    /**
     * @brief Solve the linear system Ax = b
     * @param rhs Right-hand side array b
     * @param solution Solution array x (output)
     * @param params Solver parameters
     * @throws std::invalid_argument if dimensions don't match
     */
    void solve(const utils::Array2D& rhs, utils::Array2D& solution,
              const SolverParams& params);

    /**
     * @brief Get statistics from last solve
     * @return Solver statistics
     */
    const SolverStats& get_stats() const { return stats_; }

    /**
     * @brief Get solver name
     * @return String identifier for the solver
     */
    std::string get_name() const;

    /**
     * @brief Get solver type
     * @return Solver type enum
     */
    SolverType get_type() const {
        return use_preconditioner_ ? SolverType::PreconditionedCG : SolverType::ConjugateGradient;
    }

    /**
     * @brief Reset solver state (clears statistics)
     */
    void reset();

    /**
     * @brief Set restart threshold for numerical stability
     * @param max_iterations_without_progress Max iterations without residual reduction
     *
     * If the residual doesn't decrease after this many iterations,
     * the solver restarts with current solution as initial guess.
     */
    void set_restart_threshold(int max_iterations_without_progress) {
        restart_threshold_ = max_iterations_without_progress;
    }

    /**
     * @brief Get number of restarts performed
     * @return Number of restarts
     */
    int get_restarts() const { return stats_.restarts; }

    /**
     * @brief Set lambda coefficient
     * @param lambda Diffusion coefficient dt/h²
     */
    void set_lambda(double lambda) { lambda_ = lambda; }

    /**
     * @brief Enable or disable preconditioner
     * @param use_preconditioner Enable Jacobi preconditioner
     */
    void set_preconditioner(bool use_preconditioner) {
        use_preconditioner_ = use_preconditioner;
    }

    /**
     * @brief Check if MPI is enabled
     * @return true if running in parallel
     */
    bool is_parallel() const { return is_parallel_; }

    /**
     * @brief Get MPI rank
     * @return Process rank (0 in serial mode)
     */
    int get_rank() const { return rank_; }

    /**
     * @brief Get number of processes
     * @return Number of processes (1 in serial mode)
     */
    int get_size() const { return size_; }

private:
    // Algorithm parameters
    bool use_preconditioner_;       ///< Use Jacobi preconditioner
    double lambda_;                  ///< Diffusion coefficient
    int restart_threshold_;          ///< Max iterations without progress

    // MPI parameters
    bool is_parallel_;               ///< MPI parallel execution
    MPI_Comm comm_;                  ///< MPI communicator
    int rank_;                       ///< Process rank
    int size_;                       ///< Number of processes
    int neighbor_rank_[4];           ///< Neighbor ranks [N, S, E, W]

    // Solver state
    SolverStats stats_;              ///< Statistics from last solve
    bool initialized_;               ///< Solver initialization flag

    // Work arrays for CG algorithm
    std::unique_ptr<utils::Array2D> r_;      ///< Residual vector
    std::unique_ptr<utils::Array2D> p_;      ///< Search direction
    std::unique_ptr<utils::Array2D> Ap_;     ///< Matrix-vector product
    std::unique_ptr<utils::Array2D> z_;      ///< Preconditioned residual (PCG)
    std::unique_ptr<utils::Array2D> temp_;   ///< Temporary array

    /**
     * @brief Initialize work arrays
     * @param rows Number of rows
     * @param cols Number of columns
     */
    void initialize_work_arrays(size_t rows, size_t cols);

    /**
     * @brief Matrix-vector multiplication: y = Ax
     * @param x Input vector
     * @param y Output vector y = Ax
     *
     * Implements 5-point stencil for (I - λ∇²):
     * y[i][j] = (1 + 4λ)x[i][j] - λ(x[i+1][j] + x[i-1][j] + x[i][j+1] + x[i][j-1])
     */
    void matrix_vector_multiply(const utils::Array2D& x, utils::Array2D& y);

    /**
     * @brief Compute dot product: result = x^T y
     * @param x First vector
     * @param y Second vector
     * @return Dot product (sum over all elements, reduced across MPI processes)
     */
    double dot_product(const utils::Array2D& x, const utils::Array2D& y);

    /**
     * @brief Compute dot product on local domain only
     * @param x First vector
     *ariam y Second vector
     * @return Local dot product (no MPI reduction)
     */
    double local_dot_product(const utils::Array2D& x, const utils::Array2D& y);

    /**
     * @brief Apply Jacobi preconditioner: z = M^{-1} r
     * @param r Residual vector
     * @param z Preconditioned residual
     *
     * For (I - λ∇²), the diagonal is (1 + 4λ),
     * so M^{-1} = 1/(1 + 4λ) * I
     */
    void apply_preconditioner(const utils::Array2D& r, utils::Array2D& z);

    /**
     * @brief Exchange ghost cells with neighbors
     * @param x Array with ghost cells
     *
     * Required for parallel matrix-vector multiplication
     */
    void exchange_ghost_cells(utils::Array2D& x);

    /**
     * @brief Compute L2 norm of a vector
     * @param x Input vector
     * @return L2 norm (reduced across MPI processes)
     */
    double norm(const utils::Array2D& x);

    /**
     * @brief Compute local L2 norm (no MPI reduction)
     * @param x Input vector
     * @return Local L2 norm
     */
    double local_norm(const utils::Array2D& x);

    /**
     * @brief Print iteration progress
     * @param iteration Current iteration
     * @param residual Current residual
     */
    void print_progress(int iteration, double residual) const;
};

#endif // CONJUGATE_GRADIENT_SOLVER_HPP
