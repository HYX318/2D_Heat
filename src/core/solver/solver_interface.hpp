/**
 * @file solver_interface.hpp
 * @brief Abstract base class interface for iterative solvers
 */

#ifndef SOLVER_INTERFACE_HPP
#define SOLVER_INTERFACE_HPP

#include <cstddef>
#include <string>
#include <stdexcept>
#include "../../utils/array2d.hpp"

/**
 * @class ConvergenceFailureException
 * @brief Exception thrown when solver fails to converge
 */
class ConvergenceFailureException : public std::runtime_error {
public:
    explicit ConvergenceFailureException(const std::string& message)
        : std::runtime_error(message) {}
};

/**
 * @enum SolverType
 * @brief Enumeration of supported solver types
 */
enum class SolverType {
    Jacobi,
    GaussSeidel,
    SOR,
    ConjugateGradient,
    PreconditionedCG,
    Multigrid
};

/**
 * @struct SolverParams
 * @brief Parameters for solver configuration
 */
struct SolverParams {
    double tolerance;
    size_t max_iterations;
    size_t residual_check_interval;
    bool compute_final_residual;
    double lambda;
    bool verbose;
    double omega;

    SolverParams()
        : tolerance(1e-6)
        , max_iterations(10000)
        , residual_check_interval(10)
        , compute_final_residual(true)
        , lambda(0.25)
        , verbose(false)
        , omega(1.0)
    {}

    SolverParams(double tol, size_t max_iter, size_t check_interval = 10,
                bool compute_residual = true, double lam = 0.25, bool verb = false, double w = 1.0)
        : tolerance(tol)
        , max_iterations(max_iter)
        , residual_check_interval(check_interval)
        , compute_final_residual(compute_residual)
        , lambda(lam)
        , verbose(verb)
        , omega(w)
    {}
};

/**
 * @struct SolverStats
 * @brief Statistics collected during solver execution
 */
struct SolverStats {
    bool converged;
    size_t iterations;
    double final_residual;
    double initial_residual;
    double solve_time;
    double reduction_factor;
    double residual;            ///< Current residual (for CG)
    double convergence_rate;    ///< Convergence rate (for CG)
    int restarts;               ///< Number of restarts (for CG)

    SolverStats()
        : converged(false)
        , iterations(0)
        , final_residual(0.0)
        , initial_residual(0.0)
        , solve_time(0.0)
        , reduction_factor(1.0)
        , residual(0.0)
        , convergence_rate(0.0)
        , restarts(0)
    {}

    void reset() {
        converged = false;
        iterations = 0;
        final_residual = 0.0;
        initial_residual = 0.0;
        solve_time = 0.0;
        reduction_factor = 1.0;
    }
};

/**
 * @class ISolver
 * @brief Abstract base class for iterative solvers
 */
class ISolver {
public:
    virtual ~ISolver() = default;

    virtual void solve(const utils::Array2D& rhs, utils::Array2D& solution,
                      const SolverParams& params) = 0;

    virtual SolverStats get_stats() const = 0;
    virtual std::string get_name() const = 0;
    virtual SolverType get_type() const = 0;
    virtual void reset() = 0;

protected:
    ISolver() = default;
};

#endif // SOLVER_INTERFACE_HPP
