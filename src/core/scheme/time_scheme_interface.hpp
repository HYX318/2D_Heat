/**
 * @file time_scheme_interface.hpp
 * @brief Abstract base class interface for time integration schemes
 */

#ifndef TIME_SCHEME_INTERFACE_HPP
#define TIME_SCHEME_INTERFACE_HPP

#include <cstddef>
#include <string>
#include "../../mesh/mesh2d.hpp"
#include "../solver/solver_interface.hpp"

/**
 * @enum TimeSchemeType
 * @brief Enumeration of supported time integration schemes
 */
enum class TimeSchemeType {
    ImplicitEuler,
    CrankNicolson,
    ExplicitEuler,
    RungeKutta4
};

/**
 * @struct TimeParams
 * @brief Parameters for time integration
 */
struct TimeParams {
    double dt;                ///< Time step size
    double dx;                ///< Grid spacing in x-direction
    double dy;                ///< Grid spacing in y-direction
    double alpha;             ///< Diffusion coefficient
    double t_start;           ///< Start time
    double t_end;             ///< End time
    size_t nt;                ///< Number of time steps

    TimeParams()
        : dt(0.0)
        , dx(0.0)
        , dy(0.0)
        , alpha(1.0)
        , t_start(0.0)
        , t_end(1.0)
        , nt(100)
    {}

    TimeParams(double time_step, double spacing_x, double spacing_y,
               double diffusion_coeff, double start_time, double end_time,
               size_t num_steps)
        : dt(time_step)
        , dx(spacing_x)
        , dy(spacing_y)
        , alpha(diffusion_coeff)
        , t_start(start_time)
        , t_end(end_time)
        , nt(num_steps)
    {}
};

/**
 * @struct TimeStats
 * @brief Statistics collected during time integration
 */
struct TimeStats {
    double t_current;        ///< Current time
    size_t step;              ///< Current time step
    double total_time;        ///< Total wall-clock time
    double avg_step_time;     ///< Average time per step

    TimeStats()
        : t_current(0.0)
        , step(0)
        , total_time(0.0)
        , avg_step_time(0.0)
    {}

    void reset() {
        t_current = 0.0;
        step = 0;
        total_time = 0.0;
        avg_step_time = 0.0;
    }
};

/**
 * @class ITimeScheme
 * @brief Abstract base class for time integration schemes
 *
 * This class defines the interface for time integration schemes used in
 * solving the heat equation and similar parabolic PDEs.
 *
 * All time schemes must implement:
 * - step(): Perform one time step
 * - get_stats(): Return statistics
 * - get_name(): Return scheme name
 * - get_type(): Return scheme type
 */
class ITimeScheme {
public:
    virtual ~ITimeScheme() = default;

    /**
     * @brief Perform one time step
     * @param current Current solution
     * @param next Next time step solution (output)
     * @param params Time integration parameters
     * @param solver_params Solver parameters for implicit schemes
     * @throws std::runtime_error if time step fails
     */
    virtual void step(const Mesh2D& current, Mesh2D& next,
                      const TimeParams& params,
                      const SolverParams& solver_params) = 0;

    /**
     * @brief Get statistics from time integration
     * @return TimeStats structure
     */
    virtual TimeStats get_stats() const = 0;

    /**
     * @brief Get scheme name
     * @return String name of the scheme
     */
    virtual std::string get_name() const = 0;

    /**
     * @brief Get scheme type
     * @return TimeSchemeType enum value
     */
    virtual TimeSchemeType get_type() const = 0;

    /**
     * @brief Reset the scheme to initial state
     */
    virtual void reset() = 0;

protected:
    ITimeScheme() = default;
};

#endif // TIME_SCHEME_INTERFACE_HPP
