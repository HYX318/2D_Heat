/**
 * @file interfaces_compat.hpp
 * @brief Legacy compatibility layer for Interfaces (MPI communication) functions
 *
 * This header provides backward-compatible wrappers for legacy Interfaces functions
 * that internally use the new architecture classes.
 *
 * All functions are marked with [[deprecated]] to encourage migration to the new API.
 *
 * @warning This is a legacy compatibility layer and will be removed in future versions.
 * Please migrate to the new API using GhostCellExchange and CartesianTopology.
 *
 * @see MIGRATION_GUIDE.md for migration instructions
 */

#ifndef INTERFACES_COMPAT_HPP
#define INTERFACES_COMPAT_HPP

#include <mpi.h>

// Include new architecture headers
#include "../src/mpi/ghost_cell_exchange.hpp"
#include "../src/mpi/cartesian_topology.hpp"
#include "../src/core/solver/jacobi_solver.hpp"
#include "../src/utils/array2d.hpp"

#include "../Consts.h"

// Legacy mapping for parallel Jacobi
static thread_local GhostCellExchange* legacy_ghost_exchange = nullptr;
static thread_local utils::Array2D* legacy_x0 = nullptr;
static thread_local utils::Array2D* legacy_b = nullptr;
static thread_local utils::Array2D* legacy_x = nullptr;

/**
 * @brief Exchange ghost cells with neighboring processes
 *
 * This function is called outside the iteration loop for explicit schemes.
 *
 * @param Sol Solution array with ghost cells
 * @param NeighbourRank Array of neighbor process ranks
 * @param SBD_COMM Cartesian MPI communicator
 * @param colType MPI datatype for column communication
 * @param myRank Current process rank
 * @param Nx Local grid size in x-direction
 * @param Ny Local grid size in y-direction
 *
 * @deprecated Use GhostCellExchange::exchange() with Mesh2D.
 */
[[deprecated("Use GhostCellExchange::exchange(). See MIGRATION_GUIDE()")]]
void Interfaces(double** Sol, int NeighbourRank[], MPI_Comm SBD_COMM,
               MPI_Datatype colType, int myRank, int Nx, int Ny);

/**
 * @brief Parallel Jacobi iteration solver with MPI communication
 *
 * Performs Jacobi iteration with ghost cell exchange at each iteration.
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
 * @param NeighbourRank Array of neighbor process ranks
 * @param SBD_COMM Cartesian MPI communicator
 * @param colType MPI datatype for column communication
 * @param myRank Current process rank
 *
 * @deprecated Use JacobiSolver<true> with GhostCellExchange.
 */
[[deprecated("Use JacobiSolver<true>. See MIGRATION_GUIDE()")]]
void MPIJacobi(double** x, double** x0, double** b, double& Residu, double Tol,
               int& iConv, int Nx, int Ny, double lambda,
               int NeighbourRank[], MPI_Comm SBD_COMM, MPI_Datatype colType,
               int myRank);

#endif // INTERFACES_COMPAT_HPP
