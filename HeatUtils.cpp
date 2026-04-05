#include <mpi.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include "HeatUtils.h"
#include "Interfaces.h"

/**
 * Read simulation parameters from Param.in file
 * @param Nx: Number of grid points in x-direction per sub-domain
 * @param Ny: Number of grid points in y-direction per sub-domain
 * @param Nt: Number of time steps
 * @param StabP: Stability parameter for time step calculation
 */
void ReadParam(int &Nx, int &Ny, int &Nt, double &StabP){
  std::ifstream myfile;
  myfile.open("Param.in");
  if(!myfile) {
    std::cerr << "Param.in file could not be opened" << std::endl;
  }
  myfile >> Nx >> Ny ;
  myfile >> Nt;
  myfile >> StabP;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////   2D Contiguous Array Allocation  ///////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Allocate a 2D array with contiguous memory layout
 * This is important for efficient MPI communication
 *
 * @param Tab: Pointer to 2D array
 * @param rowLength: Number of columns
 * @param colLength: Number of rows
 */
void Contiguous2D(double  **Tab, int rowLength, int colLength){
  Tab[0] = new double[colLength*rowLength];
  for (int i=1; i < colLength; i++)
    Tab[i] = Tab[0] + i * rowLength ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////     Write to File    ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Export solution to file
 * Format: I J x y solution_value
 *
 * @param Sol: 2D solution array
 * @param I: Global x-indices
 * @param J: Global y-indices
 * @param x: x-coordinates
 * @param y: y-coordinates
 * @param Nx: Local grid size in x-direction
 * @param Ny: Local grid size in y-direction
 * @param SolFile: Output filename
 */
void
Export(double **Sol, int I[], int J[], double x[], double y[], int Nx, int Ny, std::string SolFile){
  std::ofstream outfile;
  outfile.open(SolFile);
  for (int j = 1; j < Ny+1; j++)
  {
     for (int i = 1; i < Nx+1; i++){
       outfile << I[i] << " " << J[j] << " " << x[i]<< " " << y[j] << " " << std::setprecision(15) <<  Sol[j][i] << std::endl;
     }
  }
  outfile.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////          Initialization        ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Initialize solution with initial condition (t=0)
 *
 * @param Sol: Solution array to initialize
 * @param x: x-coordinates
 * @param y: y-coordinates
 * @param Nx: Local grid size in x-direction
 * @param Ny: Local grid size in y-direction
 */
void Init(double **Sol, double x[], double y[], int Nx, int Ny){
  double t0 = 0.;
  ExactSol(Sol, x, y, t0, Nx, Ny);
}

/**
 * Analytical solution at a single point
 * Solves: u_t = 5*(u_xx + u_yy) on unit square
 * With initial condition: u(x,y,0) = sin(pi*x)*sin(2*pi*y)
 * Boundary conditions: u = 0 on all boundaries
 *
 * @param x: x-coordinate
 * @param y: y-coordinate
 * @param t: time
 * @return Solution value
 */
double Analytic(double x, double y, double t){
  double pi = atan(1)*4;
  double k0 = -5*pi*pi;
  double Val = sin(pi*x)*sin(2*pi*y)*exp(t*k0);
  return Val;
}

/**
 * Compute exact solution on the entire grid
 *
 * @param Sol: Solution array
 * @param x: x-coordinates
 * @param y: y-coordinates
 * @param t: time
 * @param Nx: Local grid size in x-direction
 * @param Ny: Local grid size in y-direction
 */
void ExactSol(double **Sol, double x[], double y[], double t, int Nx, int Ny){
  for (int j=0; j<Ny+2; j++)
    for (int i=0; i<Nx+2; i++){
      Sol[j][i] = Analytic(x[i], y[j], t);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////   Error Calculation  /////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Compute L2 norm of error array
 *
 * @param Err: Error array
 * @param Nx: Local grid size in x-direction
 * @param Ny: Local grid size in y-direction
 * @param h: Grid spacing
 * @return L2 norm
 */
double TwoNorm(double **Err, int Nx, int Ny, double h){
  double e_2;
  e_2 = 0.;
  for (int j=1; j<Ny+1; j++)
      for (int i=1; i<Nx+1; i++)
        e_2 = e_2 + pow(Err[j][i], 2);
  return sqrt(e_2)*h;
}

/**
 * Compute infinity norm (maximum absolute value) of error array
 *
 * @param Err: Error array
 * @param Nx: Local grid size in x-direction
 * @param Ny: Local grid size in y-direction
 * @return Infinity norm
 */
double InftyNorm(double **Err, int Nx, int Ny){
  double e_infty;
  e_infty = 0.;
  for (int j=1; j<Ny+1; j++)
      for (int i=1; i<Nx+1; i++)
        e_infty = fmax(e_infty , fabs(Err[j][i]));
  return e_infty;
}

/**
 * Compute error between numerical and reference solutions
 *
 * @param Sol: Numerical solution
 * @param Ref: Reference/exact solution
 * @param Err: Error array (output)
 * @param Nx: Local grid size in x-direction
 * @param Ny: Local grid size in y-direction
 */
void Error(double **Sol, double **Ref, double **Err, int Nx, int Ny){
  for (int j=0; j<Ny+2; j++)
    for (int i=0; i<Nx+2; i++)
      Err[j][i] = Sol[j][i] - Ref[j][i];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////  Implicit Scheme Solver  ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Copy contents of one 2D array to another
 *
 * @param Sol: Source array
 * @param Sol0: Destination array
 * @param Nx: Local grid size in x-direction
 * @param Ny: Local grid size in y-direction
 */
void Copy(double **Sol, double **Sol0, int Nx, int Ny){
  for (int j=0; j<Ny+2; j++)
    for (int i=0; i<Nx+2; i++)
      Sol0[j][i] = Sol[j][i];
}

/**
 * Jacobi iteration solver for implicit Euler scheme
 * Solves the linear system: (I - lambda*L) * x = b
 * where L is the discrete Laplacian operator
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
 */
void Jacobi(double **x, double **x0, double **b, double &Residu, double Tol,
            int &iConv, int Nx, int Ny, double lambda){
  double Coeff;
  Coeff = double(1)/(1+4*lambda);
  double Res;
  Residu = 1.0;
  iConv = 0;
  while (Residu > Tol ){
    iConv = iConv + 1;
    Res = 0.0;
    for (int j=1; j<Ny+1; j++)
      for (int i=1; i<Nx+1; i++){
        // Jacobi update formula for 2D heat equation
        x[j][i] = Coeff * ( lambda * (x0[j+1][i] + x0[j][i+1]
                                   + x0[j-1][i] + x0[j][i-1]) + b[j][i] );
        Res = fmax(abs(x[j][i] - x0[j][i]), Res);
      }
    Residu = Res;
    Copy(x, x0, Nx, Ny);
  }
}
