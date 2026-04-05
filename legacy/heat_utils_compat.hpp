/**
 * @file heat_utils_compat.hpp
 * @brief Legacy compatibility layer for HeatUtils functions
 *
 * This header provides backward-compatible wrappers for legacy HeatUtils functions
 * that internally use the new architecture classes.
 *
 * All functions are marked with [[deprecated]] to encourage migration to the new API.
 *
 * @warning This is a legacy compatibility layer and will be removed in future versions.
 * Please migrate to the new API using Mesh2D, Array2D, and the new solver interfaces.
 *
 * @see MIGRATION_GUIDE.md for migration instructions
 */

#ifndef HEAT_UTILS_COMPAT_HPP
#define HEAT_UTILS_COMPAT_HPP

#include <mpi.h>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

// Include new architecture headers
#include "../src/utils/array2d.hpp"
#include "../src/mesh/mesh2d.hpp"
#include "../src/core/solver/jacobi_solver.hpp"

// Global objects for legacy compatibility (thread-local for thread safety)
struct LegacyCompatGlobals {
    static thread_local utils::Array2D* current_sol;
    static thread_local utils::Array2D* current_sol0;
    static thread_local utils::Array2D* current_ref;
    static thread_local utils::Array2D* current_err;
    static thread_local utils::Array2D* current_x0;
    static thread_local utils::Array2D* current_b;
};

// Static initialization
thread_local utils::Array2D* LegacyCompatGlobals::current_sol = nullptr;
thread_local utils::Array2D* LegacyCompatGlobals::current_sol0 = nullptr;
thread_local utils::Array2D* LegacyCompatGlobals::current_ref = nullptr;
thread_local utils::Array2D* LegacyCompatGlobals::current_err = nullptr;
thread_local utils::Array2D* LegacyCompatGlobals::current_x0 = nullptr;
thread_local utils::Array2D* LegacyCompatGlobals::current_b = nullptr;

/**
 * @brief Read simulation parameters from Param.in file
 * @param Nx Number of grid points in x-direction per sub-domain (output)
 * @param Ny Number of grid points in y-direction per sub-domain (output)
 * @param Nt Number of time steps (output)
 * @param StabP Stability parameter for time step calculation (output)
 *
 * @deprecated Use the new configuration system instead. See migration guide.
 */
[[deprecated("Use new configuration system. See MIGRATION_GUIDE.md")]]
void ReadParam(int& Nx, int& Ny, int& Nt, double& StabP);

/**
 * @brief Allocate a 2D array with contiguous memory layout
 * @param Tab Pointer to 2D array
 * @param rowLength Number of columns
 * @param colLength Number of rows
 *
 * @deprecated Use utils::Array2D or Mesh2D instead. See migration guide.
 *
 * Note: This wrapper allocates memory but sets up a mapping to Array2D
 */
[[deprecated("Use utils::Array2D or Mesh2D. See MIGRATION_GUIDE.md")]]
void Contiguous2D(double** Tab, int rowLength, int colLength);

/**
 * @brief Initialize solution with initial condition (t=0)
 * @param Sol Solution array to initialize
 * @param x x-coordinates
 * @param y y-coordinates
 * @param Nx Local grid size in x-direction
 * @param Ny Local grid size in y-direction
 *
 * @deprecated Use Mesh2D::apply_bc() with initial condition function.
 */
[[deprecated("Use Mesh2D::apply_bc(). See MIGRATION_GUIDE.md")]]
void Init(double** Sol, double x[], double y[], int Nx, int Ny);

/**
 * @brief Compute exact solution on the entire grid
 * @param Sol Solution array
 * @param x x-coordinates
 * @param y y-coordinates
 * @param t time
 * @param Nx Local grid size in x-direction
 * @param Ny Local grid size in y-direction
 *
 * @deprecated Use Mesh2D with analytic solution callback.
 */
[[deprecated("Use Mesh2D with analytic callback. See MIGRATION_GUIDE.md")]]
void ExactSol(double** Sol, double x[], double y[], double t, int Nx, int Ny);

/**
 * @brief Analytical solution at a single point
 * Solves: u_t = 5*(u_xx + u_yy) on unit square
 * With initial condition: u(x,y,0) = sin(pi*x)*sin(2*pi*y)
 * Boundary conditions: u = 0 on all boundaries
 *
 * @param x x-coordinate
 * @param y y-coordinate
 * @param t time
 * @return Solution value
 *
 * @deprecated Use the new analytic solution function in Mesh2D.
 */
[[deprecated("Use new analytic solution. See MIGRATION_GUIDE.md")]]
double Analytic(double x, double y, double t);

/**
 * @brief Compute error between numerical and reference solutions
 * @param Sol Numerical solution
 * @param Ref Reference/exact solution
 * @param Err Error array (output)
 * @param Nx Local grid size in x-direction
 * @param Ny Local grid size in y-direction
 *
 * @deprecated Use utils::Array2D::operator-= and norm functions.
 */
[[deprecated("Use Array2D arithmetic. See MIGRATION_GUIDE.md")]]
void Error(double** Sol, double** Ref, double** Err, int Nx, int Ny);

/**
 * @brief Compute L2 norm of error array
 * @param Err Error array
 * @param Nx Local grid size in x-direction
 * @param Ny Local grid size in y-direction
 * @param h Grid spacing
 * @return L2 norm
 *
 * @deprecated Use utils::Array2D::l2_norm().
 */
[[deprecated("Use Array2D::l2_norm(). See MIGRATION_GUIDE.md")]]
double TwoNorm(double** Err, int Nx, int Ny, double h);

/**
 * @brief Compute infinity norm (maximum absolute value) of error array
 * @param Err Error array
 * @param Nx Local grid size in x-direction
 * @param Ny Local grid size in y-direction
 * @return Infinity norm
 *
 * @deprecated Use utils::Array2D::linfty_norm().
 */
[[deprecated("Use Array2D::linfty_norm(). See MIGRATION_GUIDE.md")]]
double InftyNorm(double** Err, int Nx, int Ny);

/**
 * @brief Copy array contents
 * @param Sol Source array
 * @param Sol0 Destination array
 * @param Nx Local grid size in x-direction
 * @param Ny Local grid size in y-direction
 *
 * @deprecated Use utils::Array2D::copy_from().
 */
[[deprecated("Use Array2D::copy_from(). See MIGRATION_GUIDE.md")]]
void Copy(double** Sol, double** Sol0, int Nx, int Ny);

/**
 * @brief Jacobi iteration solver for implicit Euler scheme (serial version)
 * Solves the linear system: (I - lambda*L) * x = b
 * where L is the discrete Laplacian operator
 *
 * @param x Solution array (output)
 * @param x0 Initial guess (updated each iteration)
 * @param b Right-hand side array
 * @param Residu Residual after convergence (output)
 * @param Tol Convergence tolerance
 * @param iConv Number of iterations (output)
 * @param Nx Local grid size in x-direction
 * @param Ny Local grid size in y-direction
 * @param lambda Coefficient = Dt/(h^2)
 *
 * @deprecated Use JacobiSolver<false> for serial or JacobiSolver<true> for parallel.
 */
[[deprecated("Use JacobiSolver. See MIGRATION_GUIDE.md")]]
void Jacobi(double** x, double** x0, double** b, double& Residu, double Tol,
            int& iConv, int Nx, int Ny, double lambda);

/**
 * @brief Export solution to file
 * Format: I J x y solution_value
 *
 * @param Sol 2D solution array
 * @param I Global x-indices
 * @param J Global y-indices
 * @param x x-coordinates
 * @param y y-coordinates
 * @param Nx Local grid size in x-direction
 * @param Ny Local grid size in y-direction
 * @param SolFile Output filename
 *
 * @deprecated Use the new I/O utilities in Mesh2D.
 */
[[deprecated("Use new I/O utilities. See MIGRATION_GUIDE.md")]]
void Export(double** Sol, int I[], int J[], double x[], double y[],
            int Nx, int Ny, std::string SolFile);

#endif // HEAT_UTILS_COMPAT_HPP
