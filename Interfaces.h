void Interfaces(double **Sol, int NeighbourRank[], MPI_Comm SBD_COMM, MPI_Datatype colType, int myRank,
                int Nx, int Ny);

void MPIJacobi(double **x, double **x0, double **b, double &Residu, double Tol, int &iConv, int Nx, int Ny, double lambda, 
            int NeighbourRank[], MPI_Comm SBD_COMM, MPI_Datatype colType, int myRank);

            