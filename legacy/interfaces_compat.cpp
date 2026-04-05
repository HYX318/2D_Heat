/**
 * @file interfaces_compat.cpp
 * @brief Implementation of legacy compatibility layer for Interfaces functions
 */

#include "interfaces_compat.hpp"
#include <iostream>
#include <stdexcept>

/**
 * @brief Exchange ghost cells with neighboring processes
 *
 * This is a compatibility wrapper that maintains the same behavior as the original
 * Interfaces function but could be enhanced to use the new GhostCellExchange class.
 */
void Interfaces(double** Sol, int NeighbourRank[], MPI_Comm SBD_COMM,
               MPI_Datatype colType, int myRank, int Nx, int Ny) {
    // For now, we use the original MPI_Sendrecv pattern to maintain exact compatibility
    // In a future version, this could optionally use GhostCellExchange

    // Exchange rows (North-South)
    // Send bottom row to South, receive from North
    MPI_Sendrecv(&Sol[1][0], 1, colType, NeighbourRank[South], 0,
                 &Sol[Ny + 1][0], 1, colType, NeighbourRank[North], 0,
                 SBD_COMM, MPI_STATUS_IGNORE);

    // Send top row to North, receive from South
    MPI_Sendrecv(&Sol[Ny][0], 1, colType, NeighbourRank[North], 1,
                 &Sol[0][0], 1, colType, NeighbourRank[South], 1,
                 SBD_COMM, MPI_STATUS_IGNORE);

    // Exchange columns (East-West)
    // Send left column to West, receive from East
    MPI_Sendrecv(&Sol[0][1], Nx, MPI_DOUBLE, NeighbourRank[West], 2,
                 &Sol[0][Nx + 1], Nx, MPI_DOUBLE, NeighbourRank[East], 2,
                 SBD_COMM, MPI_STATUS_IGNORE);

    // Send right column to East, receive from West
    MPI_Sendrecv(&Sol[0][Nx], Nx, MPI_DOUBLE, NeighbourRank[East], 3,
                 &Sol[0][1], Nx, MPI_DOUBLE, NeighbourRank[West], 3,
                 SBD_COMM, MPI_STATUS_IGNORE);

    std::cout << "WARNING [legacy]: Interfaces is deprecated. "
              << "Use GhostCellExchange::exchange(). See MIGRATION_GUIDE.md" << std::endl;
}

/**
 * @brief Parallel Jacobi iteration solver with MPI communication
 *
 * This is a compatibility wrapper that maintains the same behavior as the original
 * MPIJacobi function.
 */
void MPIJacobi(double** x, double** x0, double** b, double& Residu, double Tol,
               int& iConv, int Nx, int Ny, double lambda,
               int NeighbourRank[], MPI_Comm SBD_COMM, MPI_Datatype colType,
               int myRank) {
    // Coefficient for Jac Jacobi update
    double Coeff = 1.0 / (1.0 + 4.0 * lambda);
    double Res;
    Residu = 1.0;
    iConv = 0;

    // Jacobi iteration loop
    while (Residu > Tol) {
        iConv = iConv + 1;
        Res = 0.0;

        // Jacobi update on interior points
        for (int j = 1; j < Ny + 1; j++) {
            for (int i = 1; i < Nx + 1; i++) {
                x[j][i] = Coeff * (lambda * (x0[j + 1][i] + x0[j][i + 1] +
                                              x0[j - 1][i] + x0[j][i - 1]) + b[j][i]);
                Res = std::max(std::abs(x[j][i] - x0[j][i]), Res);
            }
        }

        // Exchange ghost cells with neighboring processes
        // Exchange rows (North-South)
        MPI_Sendrecv(&x[1][0], 1, colType, NeighbourRank[South], 0,
                     &x[Ny + 1][0], 1, colType, NeighbourRank[North], 0,
                     SBD_COMM, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&x[Ny][0], 1, colType, NeighbourRank[North], 1,
                     &x[0][0], 1, colType, NeighbourRank[South], 1,
                     SBD_COMM, MPI_STATUS_IGNORE);

        // Exchange columns (East-West)
        MPI_Sendrecv(&x[0][1], Nx, MPI_DOUBLE, NeighbourRank[West], 2,
                     &x[0][Nx + 1], Nx, MPI_DOUBLE, NeighbourRank[East], 2,
                     SBD_COMM, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&x[0][Nx], Nx, MPI_DOUBLE, NeighbourRank[East], 3,
                     &x[0][1], Nx, MPI_DOUBLE, NeighbourRank[West], 3,
                     SBD_COMM, MPI_STATUS_IGNORE);

        // Copy new solution to old solution for next iteration
        for (int j = 0; j < Ny + 2; j++) {
            for (int i = 0; i < Nx + 2; i++) {
                x0[j][i] = x[j][i];
            }
        }
        Residu = Res;

        // Global reduction: compute maximum residual across all processes
        double globalResidu;
        MPI_Allreduce(&Res, &globalResidu, 1, MPI_DOUBLE, MPI_MAX, SBD_COMM);
        Residu = globalResidu;
    }

    std::cout << "WARNING [legacy]: MPIJacobi is deprecated. "
              << "Use JacobiSolver<true>. See MIGRATION_GUIDE.md" << std::endl;
}
