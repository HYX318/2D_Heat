#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
//#include <unistd.h>
// basic file operations
#include <iostream>
#include <iomanip>
#include <fstream>
// string
#include <string>
#include <algorithm>
#include "HeatUtils.h"
#include "Consts.h"
#include "Interfaces.h"

/**
 * 2D Heat Equation Solver using MPI Parallel Computing
 * Implements implicit Euler scheme with Jacobi iteration
 *
 * Domain decomposition: 3x3 Cartesian grid (9 MPI processes required)
 */
int main(int argc,char **argv)
{
  // MPI initialization
  MPI_Init(&argc, &argv);
  int myRank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // Process grid configuration (3x3 decomposition)
  const int Npx = 3, Npy = 3;

  // Check if number of processes matches decomposition
  if (numProcs != Npx * Npy) {
      if (myRank == 0) {
          std::cerr << "Error: Required " << Npx * Npy << " processes, but got " << numProcs << " processes." << std::endl;
      }
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  // Parameter definition
  int Nx, Ny, Nnt;
  double StabP, Dt;

  // Read parameters only on root process
  if (myRank == 0) {
    ReadParam(Nx, Ny, Nnt, StabP);
  }

  // Broadcast parameters from root to all processes
  MPI_Bcast(&Nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Nnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&StabP, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Mesh Definition
  int Nx_tot, Ny_tot;
  Nx_tot = Nx; Ny_tot = Ny;

  // Compute grid spacing and time step
  double hx = 1.0/(Nx_tot+1);
  double hy = 1.0/(Ny_tot+1);
  double h = hx;
  Dt = StabP*(h*h);

  // Create Cartesian Communicator for domain decomposition
  MPI_Comm SBD_COMM;
  const int nDims = 2;
  int dims[nDims] = {Npx, Npy};
  int periods[nDims] = {0, 0}; // Non-periodic boundaries
  MPI_Cart_create(MPI_COMM_WORLD, nDims, dims, periods, 1, &SBD_COMM);

  // Get rank in Cartesian communicator
  int cartRank;
  MPI_Comm_rank(SBD_COMM, &cartRank);

  // Get process coordinates in Cartesian grid
  int myCoord[nDims];
  MPI_Cart_coords(SBD_COMM, myRank, nDims, myCoord);

  // Get neighbor process ranks
  int NeighbourRank[2*nDims];
  // X direction (West-East)
  MPI_Cart_shift(SBD_COMM, 0, 1, &NeighbourRank[West], &NeighbourRank[East]);
  // Y direction (South-North)
  MPI_Cart_shift(SBD_COMM, 1, 1, &NeighbourRank[South], &NeighbourRank[North]);

  // Dump process information to file
  std::ofstream procFile("Proc-" + std::to_string(myRank) + ".txt");
  procFile << "My Rank: " << myRank << "\n";
  procFile << "My Coordinates: (" << myCoord[0] << ", " << myCoord[1] << ")\n";
  procFile << "Neighbours: " << "\n";
  procFile << "South (" << NeighbourRank[South] << ")\n";
  procFile << "North (" << NeighbourRank[North] << ")\n";
  procFile << "West (" << NeighbourRank[West] << ")\n";
  procFile << "East (" << NeighbourRank[East] << ")\n";
  procFile.close();

  // Compute local domain boundaries
  int iMin, iMax, jMin, jMax;
  iMin = 1; jMin = 1;
  iMax = Nx_tot; jMax = Ny_tot;

  // Get Cartesian grid information
  MPI_Cart_get(SBD_COMM, 2, dims, periods, myCoord);

  // Compute local index range for this sub-domain
  iMin = myCoord[0] * Nx + 1;
  iMax = iMin + Nx - 1;
  jMin = myCoord[1] * Ny + 1;
  jMax = jMin + Ny - 1;

  // Allocate arrays with ghost cells (colLength = Ny+2, rowLength = Nx+2)
  int colLength = Ny+2; int rowLength = Nx+2;

  // Solution arrays
  double** Sol = new double*[colLength];
  Contiguous2D(Sol, colLength, rowLength);
  double** Sol0 = new double*[colLength];
  Contiguous2D(Sol0, colLength, rowLength);

  // Reference/exact solution array
  double** Ref = new double*[colLength];
  Contiguous2D(Ref, colLength, rowLength);

  // Error array
  double** Err = new double*[colLength];
  Contiguous2D(Err, colLength, rowLength);

  // Coordinate and index arrays
  double* x = new double[Nx+2];
  double* y = new double[Ny+2];
  int* I = new int[Nx+2];
  int* J = new int[Ny+2];

  // Compute global coordinates and indices
  for (int i=0; i<Nx+2; i++) {
    I[i] = iMin-1+i;
    x[i] = I[i] * hx;
  }
  for (int j=0; j<Ny+2; j++) {
    J[j] = jMin-1+j;
    y[j] = J[j] * hy;
  }

  // Initialize solution with initial condition
  double T = 0;
  Init(Sol, x, y, Nx, Ny);

  // Arrays for implicit scheme solver
  double** x0 = new double*[colLength]; // Initial guess for iterative solver
  Contiguous2D(x0, colLength, rowLength);
  double** b = new double*[colLength]; // RHS of linear system
  Contiguous2D(b, colLength, rowLength);

  // Implicit scheme parameters
  double Residu = 1.0;
  double DtImp = Dt*100;
  double lambda = DtImp/(h*h);
  double Tol = 1.e-6;
  int iConv;

  // Create MPI datatype for column communication (vector type)
  MPI_Datatype colType;
  MPI_Type_vector(colLength, 1, rowLength, MPI_DOUBLE, &colType);
  MPI_Type_commit(&colType);

  // Time stepping loop
  for (int nStep = 0; nStep < Nnt; nStep++){
    T = T + DtImp;
    iConv = 0;

    // Set RHS of linear system
    Copy(Sol, b, Nx, Ny);

    // Set initial guess for Jacobi iteration
    Copy(Sol, x0, Nx, Ny);

    // Solve implicit scheme using parallel Jacobi iteration
    // Note: Use MPIJacobi for parallel, Jacobi for serial
    MPIJacobi(Sol, x0, b, Residu, Tol, iConv, Nx, Ny, lambda, NeighbourRank, SBD_COMM, colType, myRank);

    // Display iteration info (only from root process)
    if (myRank == 0) {
      std::cout << "Step=" << nStep << " Residu=" << Residu << " Iter=" << iConv << std::endl;
    }
  }

  // Compute exact solution at final time
  ExactSol(Ref, x, y, T, Nx, Ny);

  // Compute error between numerical and exact solutions
  Error(Sol, Ref, Err, Nx, Ny);

  // Compute local error norms
  double e2, einfty;
  e2 = pow(TwoNorm(Err, Nx, Ny, h), 2);
  einfty = InftyNorm(Err, Nx, Ny);

  // Reduce error norms across all processes
  double Global_e2, Global_einfty;
  MPI_Reduce(&e2, &Global_e2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&einfty, &Global_einfty, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  // Display results (only from root process)
  if (myRank == 0) {
    std::cout << "Numerical Approximation: [Euler Implicit]" << std::endl;
    std::cout << " ## Tol=" << Tol << " StabP=" << DtImp/(h*h) << std::endl;
    std::cout << " ** N=" << Nx << " N_tot=" << Ny_tot << std::endl;
    std::cout << " ** h=" << hx << " Dt=" << DtImp << std::endl;
    std::cout << " ** T=" << T << " T/Dt=" << T/DtImp << std::endl;
    std::cout << " ** Error L2-norm: " << sqrt(Global_e2) << std::endl;
    std::cout << " ** Error Li-norm: " << Global_einfty << std::endl;

    // Export results to files
    std::string ThisFile = "ExactImp.txt";
    Export(Ref, I, J, x, y, Nx, Ny, ThisFile);
    ThisFile = "SolImp.txt";
    Export(Sol, I, J, x, y, Nx, Ny, ThisFile);
  }

  // Free memory
  delete[] Sol[0]; delete[] Sol;
  delete[] Sol0[0]; delete[] Sol0;
  delete[] Ref[0]; delete[] Ref;
  delete[] Err[0]; delete[] Err;
  delete[] x0[0]; delete[] x0;
  delete[] b[0]; delete[] b;
  delete[] x; delete[] y; delete[] I; delete[] J;

  // Free MPI datatype
  MPI_Type_free(&colType);

  // Terminate MPI
  MPI_Finalize();

  return 0;
}
