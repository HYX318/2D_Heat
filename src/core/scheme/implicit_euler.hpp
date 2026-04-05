/**
 * @file implicit_euler.hpp
 * @brief Implicit Euler time integration scheme for the heat equation
 */

#ifndef IMPLICIT_EULER_HPP
#define IMPLICIT_EULER_HPP

#include "time_scheme_interface.hpp"
#include "../../utils/logger.hpp"

/**
 * @class ImplicitEuler
 * @brief Implicit Euler time integration scheme
 *
 * The implicit Euler scheme for the heat equation ∂u/∂t = α∇²u is:
 * (u^{n+1} - u^n) / Δt = α ∇²u^{n+1}
 *
 * Rearranging to solve for u^{n+1}:
 * u^{n+1} - αΔt ∇²u^{n+1} = u^n
 *
 * This requires solving a linear system at each time step, which is
 * unconditionally stable but requires an iterative solver.
 *
 * Features:
 * - Unconditionally stable for any Δt
 * - First-order accurate in time
 * - Requires solving linear system
 * - Supports any iterative solver through interface
 */
class ImplicitEuler : public ITimeScheme {
public:
    /**
     * @brief Constructor
     */
    ImplicitEuler();

    /**
     * @brief Destructor
     */
    virtual ~ImplicitEuler() = default;

    /**
     * @brief Perform one time step using implicit Euler
     * @param current Current solution at time t^n
     * @param next Next time step solution at time t^{n+1} (output)
     * @param params Time integration parameters
     * @param solver_params Parameters for the linear solver
     * @throws std::runtime_error if time step fails
     *
     * The method:
     * 1. Constructs the RHS: (1/Δt) * u^n
     * 2. Solves (I - αΔt ∇²) u^{n+1} = (1/Δt) * u^n
     * 3. Updates statistics
     */
    void step(const Mesh2D& current, Mesh2D& next,
              const TimeParams& params,
              const SolverParams& solver_params) override;

    /**
     * @brief Get statistics from time integration
     * @return TimeStats structure
     */
    TimeStats get_stats() const override { return stats_; }

    /**
     * @brief Get scheme name
     * @return "ImplicitEuler"
     */
    std::string get_name() const override { return "ImplicitEuler"; }

    /**
     * @brief Get scheme type
     * @return TimeSchemeType::ImplicitEuler
     */
    TimeSchemeType get_type() const override { return TimeSchemeType::ImplicitEuler; }

    /**
     * @brief Reset the scheme to initial state
     */
    void reset() override {
        stats_.reset();
    }

private:
    TimeStats stats_;  ///< Statistics for this scheme
    std::unique_ptr<ISolver> solver_;  ///< Iterative solver for linear system

    /**
     * @brief Initialize the solver with appropriate parameters
     * @param params Time integration parameters
     * @param solver_params User-provided solver parameters
     */
    void initialize_solver(const TimeParams& params,
                           const SolverParams& solver_params);
};

#endif // IMPLICIT_EULER_HPP
