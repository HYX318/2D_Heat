
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "Interfaces.h"
#include "Consts.h"
#include "HeatUtils.h"
void MPIJacobi(double **x, double **x0, double **b, double &Residu, double Tol, int &iConv, int Nx, int Ny, double lambda, 
            int NeighbourRank[], MPI_Comm SBD_COMM, MPI_Datatype colType, int myRank){
  double Coeff;
  Coeff = double(1)/(1+4*lambda);
  double Res;
  Residu = 1.0;
  iConv = 0;
  while (Residu > Tol ){
    iConv = iConv + 1;
    Res= 0.0;
    for (int j=1; j<Ny+1; j++)
      for (int i=1; i<Nx+1; i++){
        x[j][i] = Coeff* ( lambda* (x0[j+1][i]+x0[j][i+1]
                                   +x0[j-1][i]+x0[j][i-1]) + b[j][i] );
        Res = fmax(abs(x[j][i] - x0[j][i]),Res);
    }

    // Task 2: Update interface nodes.
    MPI_Sendrecv(&x[1][0], 1, colType, NeighbourRank[0], 0, &x[Ny + 1][0], 1, colType, NeighbourRank[1], 0, SBD_COMM, MPI_STATUS_IGNORE); // Left-Right
    MPI_Sendrecv(&x[Ny][0], 1, colType, NeighbourRank[1], 1, &x[0][0], 1, colType, NeighbourRank[0], 1, SBD_COMM, MPI_STATUS_IGNORE);   // Right-Left
    MPI_Sendrecv(&x[0][1], Nx, MPI_DOUBLE, NeighbourRank[2], 2, &x[Ny + 1][1], Nx, MPI_DOUBLE, NeighbourRank[3], 2, SBD_COMM, MPI_STATUS_IGNORE); // Top-Bottom
    MPI_Sendrecv(&x[Ny][1], Nx, MPI_DOUBLE, NeighbourRank[3], 3, &x[0][1], Nx, MPI_DOUBLE, NeighbourRank[2], 3, SBD_COMM, MPI_STATUS_IGNORE);   // Bottom-Top

    // Copy New x^{k+1} into x^k
    Copy(x,x0,Nx,Ny);
    Residu = Res;
    // Task 3: Global Residu: Reduce Residu and distribute the result back to all processes
    double globalResidu;
    MPI_Allreduce(&Res, &globalResidu, 1, MPI_DOUBLE, MPI_MAX, SBD_COMM);
    Residu = globalResidu;
  }
}
