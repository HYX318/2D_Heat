#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "Interfaces.h"
#include "Consts.h"
#include "HeatUtils.h"

/**
 * Parallel Jacobi iteration solver with MPI communication
 * This is the core parallel solver that exchanges ghost cells
 * with neighboring processes at each iteration
 *
 * @param x: Solution array (output)
 * @param x0: Initial guess (updated each iteration)
 * @param b: Right-hand side array
 * @param Residu: Residual after convergence (output)
 * @param Tol: Convergence tolerance
 * @param iConv: Number of iterations (output)
 * @param Nx: Local grid size in x-direction
 * @param Ny: Local grid size in y-direction
 * @param lambda: Coefficient = Dt/(h^2)
 * @param NeighbourRank: Array of neighbor process ranks
 * @param SBD_COMM: Cartesian MPI communicator
 * @param colType: MPI datatype for columnar communication
 * @param myRank: Current process rank
 */
void MPIJacobi(double **x, double **x0, double **b, double &Residu, double Tol, int &iConv, int Nx, int Ny, double lambda,
            int NeighbourRank[], MPI_Comm SBD_COMM, MPI_Datatype colType, int myRank){
  // Coefficient for Jacobi update
  double Coeff;
  Coeff = double(1)/(1+4*lambda);
  double Res;
  Residu = 1.0;
  iConv = 0;

  // Jacobi iteration loop
  while (Residu > Tol ){
    iConv = iConv + 1;
    Res = 0.0;

    // Jacobi update on interior points
    for (int j=1; j<Ny+1; j++)
      for (int i=1; i<Nx+1; i++){
        x[j][i] = Coeff * ( lambda * (x0[j+1][i] + x0[j][i+1]
                                   + x0[j-1][i] + x0[j][i-1]) + b[j][i] );
        Res = fmax(abs(x[j][i] - x0[j][i]), Res);
      }

    /**
     * Exchange ghost cells with neighboring processes
     *
     * MPI_Sendrecv combines send and receive in one call,
     * avoiding potential deadlock
     *
     * Data layout:
     * - x[0][i]    : South ghost row (bottom)
     * - x[Ny+1][i] : North ghost row (top)
     * - x[j][0]    : West ghost column (left)
     * - x[j][Nx+1] : East ghost column (right)
     */

    // Exchange rows (North-South)
    // Send bottom row to South, receive from North
    MPI_Sendrecv(&x[1][0], 1, colType, NeighbourRank[South], 0,
                 &x[Ny+1][0], 1, colType, NeighbourRank[North], 0,
                 SBD_COMM, MPI_STATUS_IGNORE);

    // Send top row to North, receive from South
    MPI_Sendrecv(&x[Ny][0], 1, colType, NeighbourRank[North], 1,
                 &x[0][0], 1, colType, NeighbourRank[South], 1,
                 SBD_COMM, MPI_STATUS_IGNORE);

    // Exchange columns (East-West)
    // Send left column to West, receive from East
    MPI_Sendrecv(&x[0][1], Nx, MPI_DOUBLE, NeighbourRank[West], 2,
                 &x[0][Nx+1], Nx, MPI_DOUBLE, NeighbourRank[East], 2,
                 SBD_COMM, MPI_STATUS_IGNORE);

    // Send right column to East, receive from West
    MPI_Sendrecv(&x[0][Nx], Nx, MPI_DOUBLE, NeighbourRank[East], 3,
                 &x[0][1], Nx, MPI_DOUBLE, NeighbourRank[West], 3,
                 SBD_COMM, MPI_STATUS_IGNORE);

    // Copy new solution to old solution for next iteration
    Copy(x, x0, Nx, Ny);
    Residu = Res;

    // Global reduction: compute maximum residual across all processes
    // This ensures convergence is checked globally, not just locally
    double globalResidu;
    MPI_Allreduce(&Res, &globalResidu, 1, MPI_DOUBLE, MPI_MAX, SBD_COMM);
    Residu = globalResidu;
  }
}

/**
 * Exchange ghost cells with neighboring processes
 * This function is called outside the iteration loop for explicit schemes
 *
 * @param Sol: Solution array with ghost cells
 * @param NeighbourRank: Array of neighbor process ranks
 * @param SBD_COMM: Cartesian MPI communicator
 * @param colType: MPI datatype for columnar communication
 * @param myRank: Current process rank
 * @param Nx: Local grid size in x-direction
 * @param Ny: Local grid size in y-direction
 */
void Interfaces(double **Sol, int NeighbourRank[], MPI_Comm SBD_COMM, MPI_Datatype colType, int myRank,
                int Nx, int Ny){
  // Exchange rows (North-South)
  MPI_Sendrecv(&Sol[1][0], 1, colType, NeighbourRank[South], 0,
               &Sol[Ny+1][0], 1, colType, NeighbourRank[North], 0,
               SBD_COMM, MPI_STATUS_IGNORE);
  MPI_Sendrecv(&Sol[Ny][0], 1, colType, NeighbourRank[North], 1,
               &Sol[0][0], 1, colType, NeighbourRank[South], 1,
               SBD_COMM, MPI_STATUS_IGNORE);

  // Exchange columns (East-West)
  MPI_Sendrecv(&Sol[0][1], Nx, MPI_DOUBLE, NeighbourRank[West], 2,
               &Sol[0][Nx+1], Nx, MPI_DOUBLE, NeighbourRank[East], 2,
               SBD_COMM, MPI_STATUS_IGNORE);
  MPI_Sendrecv(&Sol[0][Nx], Nx, MPI_DOUBLE, NeighbourRank[East], 3,
               &Sol[0][1], Nx, MPI_DOUBLE, NeighbourRank[West], 3,
               SBD_COMM, MPI_STATUS_IGNORE);
}
