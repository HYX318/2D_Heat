/**
 * HeatUtils.h - Utility functions for 2D Heat Equation Solver
 */

// Read simulation parameters from input file
void ReadParam(int &Nx, int &Ny, int &Nt, double &StabP);

// Allocate contiguous 2D array
void Contiguous2D(double **Tab, int rowLength, int colLength);

// Compute exact solution on the grid
void ExactSol(double **Sol, double x[], double y[], double t, int Nx, int Ny);

// Analytical solution at a single point
double Analytic(double x, double y, double t);

// Initialize solution with initial condition
void Init(double **Sol, double x[], double y[], int Nx, int Ny);

// Compute L2 norm of error
double TwoNorm(double **Err, int Nx, int Ny, double h);

// Compute infinity norm of error
double InftyNorm(double **Err, int Nx, int Ny);

// Compute error between two solutions
void Error(double **Sol, double **Ref, double **Err, int Nx, int Ny);

// Export solution to file
void Export(double **Sol, int I[], int J[], double x[], double y[], int Nx, int Ny, std::string SolFile);

// Copy array contents
void Copy(double **Sol, double **Sol0, int Nx, int Ny);

// Jacobi iteration solver for implicit scheme
void Jacobi(double **x, double **x0, double **b, double &Residu, double Tol,
            int &iConv, int Nx, int Ny, double lambda);
