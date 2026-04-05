/**
 * Interfaces.h - MPI communication functions for domain decomposition
 */

/**
 * Exchange ghost cell data with neighboring processes
 *
 * @param Sol: Solution array with ghost cells
 * @param NeighbourRank: Array of neighbor process ranks
 * @param SBD_COMM: Cartesian MPI communicator
 * @param colType: MPI datatype for column communication
 * @param myRank: Current process rank
 * @param Nx: Local grid size in x-direction
 * @param Ny: Local grid size in y-direction
 */
void Interfaces(double **Sol, int NeighbourRank[], MPI_Comm SBD_COMM, MPI_Datatype colType, int myRank,
                int Nx, int Ny);

/**
 * Parallel Jacobi iteration solver with MPI communication
 * Performs Jacobi iteration with ghost cell exchange at each iteration
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
 * @param colType: MPI datatype for column communication
 * @param myRank: Current process rank
 */
void MPIJacobi(double **x, double **x0, double **b, double &Residu, double Tol, int &iConv, int Nx, int Ny, double lambda,
            int NeighbourRank[], MPI_Comm SBD_COMM, MPI_Datatype colType, int myRank);
