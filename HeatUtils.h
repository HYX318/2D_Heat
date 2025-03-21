void ReadParam(int &Nx, int &Ny, int &Nt, double &StabP);

void Contiguous2D(double **Tab, int rowLength, int colLength);

void ExactSol(double **Sol, double x[], double y[], double t, int Nx, int Ny);

double Analytic(double x, double y, double t);

void Init(double **Sol, double x[], double y[], int Nx, int Ny);

double TwoNorm(double **Err, int Nx, int Ny, double h);

double InftyNorm(double **Err, int Nx, int Ny);

void Error(double **Sol, double **Ref, double **Err, int Nx, int Ny);

void Export(double **Sol, int I[], int J[], double x[], double y[], int Nx, int Ny, std::string SolFile);

void Copy(double **Sol, double **Sol0, int Nx, int Ny);

void Jacobi(double **x, double **x0, double **b, double &Residu, double Tol, 
            int &iConv, int Nx, int Ny, double lambda);