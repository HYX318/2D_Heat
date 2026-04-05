/**
 * @file crank_nicolson.hpp
 * @brief Crank-Nicolson time integration scheme for the heat equation
 */

#ifndef CRANK_NICOLSON_HPP
#define CRANK_NICOLSON_HPP

#include "time_scheme_interface.hpp"
#include "../../utils/logger.hpp"

/**
 * @class CrankNicolson
 * @brief Crank-Nicolson time integration scheme
 *
 * The Crank-Nicolson scheme for the heat equation ∂u/∂t = α∇²u is:
 * (u^{n+1} - u^n) / Δt = α/2 (∇²u^{n+1} + ∇²u^n)
 *
 * Rearranging to solve for u^{n+1}:
 * u^{n+1} - (αΔt/2) ∇²u^{n+1} = u^n + (αΔt/2) ∇²u^n
 *
 * This requires solving a linear system at each time step.
 *
 * Features:
 * - Unconditionally stable for any Δt
 * - Second-order accurate in time
 * - Requires solving linear system
 * - Better accuracy than implicit Euler with similar cost
 */
class CrankNicolson : public ITimeScheme {
public:
    /**
     * @brief Constructor
     */
    CrankNicolson();

    /**
     * @brief Destructor
     */
    virtual ~CrankNicolson() = default;

    /**
     * @brief Perform one time step using Crank-Nicolson
     * @param current Current solution at time t^n
     * @param next Next time step solution at time t^{n+1} (output)
     * @param params Time integration parameters
     * @param solver_params Parameters for the linear solver
     * @throws std::runtime_error if time step fails
     *
     * The method:
     * 1. Computes ∇²u^n
     * 2. Constructs the RHS: u^n + (αΔt/2) ∇²u^n
     * 3. Solves (I - αΔt/2 ∇²) u^{n+1} = RHS
     * 4. Updates statistics
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
     * @return "CrankNicolson"
     */
    std::string get_name() const override { return "CrankNicolson"; }

    /**
     * @brief Get scheme type
     * @return TimeSchemeType::CrankNicolson
     */
    TimeSchemeType get_type() const override { return TimeSchemeType::CrankNicolson; }

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

#endif // CRANK_NICOLSON_HPP
