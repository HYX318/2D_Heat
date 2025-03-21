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
int main(int argc,char **argv)
{
  // Param definition
  int Nx, Ny, Nt;
  const int Npx=3, Npy=3;
  double StabP, Dt;
  ReadParam(Nx,Ny,Nt,StabP);
  // Task 1: MPI environment
  MPI_Init(&argc, &argv);
  int myRank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // Task 1: Check MPI environment consistency
  if (numProcs != Npx * Npy) {
      if (myRank == 0) {
          std::cerr << "Error: Required " << Npx * Npy << " processes, but got " << numProcs << " processes." << std::endl;
      }
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  // Task 1: Broadcast params: Nx, Ny, Nt, StabP
  MPI_Bcast(&Nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Nt, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&StabP, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Mesh Definition
  int Nx_tot, Ny_tot;
  Nx_tot = Nx ; Ny_tot = Ny;

  // Task 1: Compute Nx_tot and Ny_tot
  double hx = 1.0/(Nx_tot+1); 
  double hy = 1.0/(Ny_tot+1);
  double h = hx;
  Dt = StabP*(h*h);
  
  

  // Task 2: Create Cartesian Communicator: SBD_COMM
  MPI_Comm SBD_COMM;
  const int nDims = 2;
  int dims[nDims] = {Npx, Npy};
  int periods[nDims] = {0, 0}; // Non-periodic boundaries
  MPI_Cart_create(MPI_COMM_WORLD, nDims, dims, periods, 1, &SBD_COMM);

  // Task 2: Get rank in SBD_COMM
  int cartRank;
  MPI_Comm_rank(SBD_COMM, &cartRank);

  // Get Coordinates in SBD_COMM
  //const int nDims = 2;
  int myCoord[nDims];
  MPI_Cart_coords(SBD_COMM, myRank, nDims, myCoord);

  // Get Neighbours
  int NeighbourRank[2*nDims];
  // X方向（东西）
  MPI_Cart_shift(SBD_COMM, 0, 1, &NeighbourRank[West], &NeighbourRank[East]);
  // Y方向（南北）
  MPI_Cart_shift(SBD_COMM, 1, 1, &NeighbourRank[South], &NeighbourRank[North]);


  // Dump to disk the file Proc-X.txt
  std::ofstream procFile("Proc-" + std::to_string(myRank) + ".txt");
  procFile << "My Rank: " << myRank << "\n";
  procFile << "My Coordinates: (" << myCoord[0] << ", " << myCoord[1] << ")\n";

  procFile << "Neighbours: " << "\n";
  procFile << "South (" << NeighbourRank[South] << ")\n";
  procFile << "North (" << NeighbourRank[North] << ")\n";
  procFile << "West (" << NeighbourRank[West] << ")\n";
  procFile << "East (" << NeighbourRank[East] << ")\n";

  procFile.close();

  
  // Node numbering
  int iMin, iMax, jMin, jMax;
  iMin = 1 ; jMin = 1 ;
  iMax = Nx_tot ; jMax = Ny_tot ;

  // MPI Node numbering
  //int myCoord[2];
  MPI_Cart_get(SBD_COMM, 2, dims, periods, myCoord);

  // Global index computation
  iMin = myCoord[0] * Nx + 1;
  iMax = iMin + Nx - 1;
  jMin = myCoord[1] * Ny + 1;
  jMax = jMin + Ny - 1;
  
   // Common (physical) row and column dimensions for all sub-domains
  int colLength = Ny+2; int rowLength = Nx+2;
  // 2D contiguous array allocation
  double** Sol = new double*[colLength];
  Contiguous2D(Sol,colLength,rowLength);
  double** Sol0 = new double*[colLength];
  Contiguous2D(Sol0,colLength,rowLength);
  // Exact Solution 
  double** Ref= new double*[colLength];
  Contiguous2D(Ref,colLength,rowLength);
  // Error
  double** Err= new double*[colLength];
  Contiguous2D(Err,colLength,rowLength);
  // Node position and numbering
  double* x = new double[Nx+2];
  double* y = new double[Ny+2];
  int* I = new int[Nx+2];
  int* J = new int[Ny+2];
  for (int i=0; i<Nx+2; i++) {
    I[i] = iMin-1+i;
    x[i] = I[i] * hx; 
  }
  for (int j=0; j<Ny+2; j++) {
    J[j] = jMin-1+j;
    y[j] = J[j] * hy;
    
  }

  // Init
  double T = 0;
  Init(Sol,x,y,Nx,Ny);
  //  Implicit Scheme
  double** x0=new double*[colLength]; // Guess for the itrative Solver
  Contiguous2D(x0,colLength,rowLength);
  double** b=new double*[colLength]; // Linear System RHS
  Contiguous2D(b,colLength,rowLength);

  double Residu=1.0;
  double DtImp = Dt*100;
  double lambda;
  lambda = DtImp/(h*h);
  double Tol = 1.e-6;
  int iConv;

  // Task 2: Create colType to exchange columns between neighbours
  MPI_Datatype colType;
  MPI_Type_vector(colLength, 1, rowLength, MPI_DOUBLE, &colType);
  MPI_Type_commit(&colType);
  
  for (int Nt = 0; Nt < 5; Nt++){
    T = T + DtImp;
    iConv = 0;
    // Rhs of the linear system
    Copy(Sol,b,Nx,Ny);    
    // Initial Guess for the Jacobi iterations
    Copy(Sol,x0,Nx,Ny);
    // Implicit Scheme: Jacobi Iteration
    Jacobi(Sol,x0,b,Residu,Tol,iConv,Nx,Ny,lambda);
    //MPIJacobi(Sol,x0,b,Residu,Tol,iConv,Nx,Ny,lambda,NeighbourRank,SBD_COMM,colType,myRank);
    // Task 2: Display messages only by one MPI proc
      std::cout << "it="<< Nt <<" Residu=" << Residu << " k="
        << iConv << std::endl;  
  }
  // Dump Analytic Solution and Numerical approximation
  // Task 2: Transform in Dumping Local solutions.
  ExactSol(Ref,x,y,T,Nx,Ny);
  std::string ThisFile = "ExactImp.txt";
  Export(Ref,I,J,x,y,Nx,Ny,ThisFile); 
  ThisFile = "SolImp.txt";
  Export(Sol,I,J,x,y,Nx,Ny,ThisFile);

  // Compute error norms
  Error(Sol,Ref,Err,Nx,Ny);
  double e2, einfty;

  e2 = pow(TwoNorm(Err,Nx,Ny,h),2);
  einfty = InftyNorm(Err,Nx,Ny);
  // Task 2: Reduce Error norms
  double Global_e2, Global_einfty;
  Global_e2 = e2; 
  Global_einfty = einfty;
  

  // Task 2: Display message only by one MPI Proc
    std::cout << "Numerical Approximation: [Euler Implicit]" <<std::endl;   
    std::cout << " ## Tol=" << Tol << " StabP=" << DtImp/(h*h) << std::endl;
    std::cout << " ** N=" <<Nx <<" N_tot="<<Ny_tot<< std::endl;
    std::cout << " ** h=" <<hx <<" Dt="<<DtImp<< std::endl; 
    std::cout << " ** T=" <<T <<" T/Dt="<<T/DtImp<< std::endl;  
    std::cout << " ** Error L2-norm: " << sqrt(Global_e2) << std::endl;
    std::cout << " ** Error Li-norm: " << Global_einfty << std::endl;
  
   // Free memory 
  delete[] Sol[0]; delete[] Sol;
  delete[] Sol0[0]; delete[] Sol0;

  delete[] Ref[0]; delete[] Ref;
  delete[] Err[0]; delete[] Err;
  delete[] x0[0]; delete[] x0;
  delete[] b[0]; delete[] b;
  delete[] x; delete[] y; delete[] I; delete[] J;
  
  // Task 2: Free memomry
  MPI_Type_free(&colType);

  // Task 1: Terminate MPI
  MPI_Finalize();
  
  return 0;
}
