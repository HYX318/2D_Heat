/**
 * @file heat_utils_compat.cpp
 * @brief Implementation of legacy compatibility layer for HeatUtils functions
 */

#include "heat_utils_compat.hpp"
#include <unordered_map>

// Mapping from double** to Array2D objects
static std::unordered_map<double**, utils::Array2D*> array_map;

/**
 * @brief Read simulation parameters from Param.in file
 */
void ReadParam(int& Nx, int& Ny, int& Nt, double& StabP) {
    std::ifstream myfile;
    myfile.open("Param.in");
    if (!myfile) {
        std::cerr << "WARNING [legacy]: Param.in file could not be opened" << std::endl;
        std::cerr << "WARNING [legacy]: Using default parameters Nx=50, Ny=50, Nt=100, StabP=0.25" << std::endl;
        Nx = 50;
        Ny = 50;
        Nt = 100;
        StabP = 0.25;
        return;
    }
    myfile myfile >> Nx >> Ny;
    myfile myfile >> Nt;
    myfile myfile >> StabP;
    myfile.close();

    std::cout << "WARNING [legacy]: ReadParam is deprecated. "
              << "Use new configuration system. See MIGRATION_GUIDE.md" << std::endl;
}

/**
 * @brief Allocate a 2D array with contiguous memory layout
 *
 * This wrapper creates an Array2D object and maps the double** pointer to it.
 * The legacy code can still use double** but the memory is managed by Array2D.
 */
void Contiguous2D(double** Tab, int rowLength, int colLength) {
    // Create Array2D with ghost cells (colLength = Ny+2, rowLength = Nx+2)
    utils::Array2D* array = new utils::Array2D(colLength, rowLength);
    array_map[Tab] = array;

    // Set up double** pointers to point into Array2D's contiguous memory
    // Note: Array2D uses row-major layout: data_[i * cols + j]
    // Legacy code uses Tab[j][i] where j is row, i is column
    // So we map Tab[j] to data row j
    double* base = array->data();
    for (int i = 0; i < colLength; i++) {
        Tab[i] = base + i * rowLength;
    }

    std::cout << "WARNING [legacy]: Contiguous2D is deprecated. "
              << "Use utils::Array2D or Mesh2D. See MIGRATION_GUIDE.md" << std::endl;
}

/**
 * @brief Analytical solution at a single point
 */
double Analytic(double x, double y, double t) {
    double pi = atan(1) * 4;
    double k0 = -5 * pi * pi;
    double Val = sin(pi * x) * sin(2 * pi * y) * exp(t * k0);
    return Val;
}

/**
 * @brief Compute exact solution on the entire grid
 */
void ExactSol(double** Sol, double x[], double y[], double t, int Nx, int Ny) {
    auto it = array_map.find(Sol);
    if (it == array_map.end()) {
        std::cerr << "ERROR [legacy]: ExactSol called with unmapped array" << std::endl;
        return;
    }
    utils::Array2D* array = it->second;

    for (int j = 0; j < Ny + 2; j++) {
        for (int i = 0; i < Nx + 2; i++) {
            Sol[j][i] = Analytic(x[i], y[j], t);
        }
    }
}

/**
 * @brief Initialize solution with initial condition
 */
void Init(double** Sol, double x[], double y[], int Nx, int Ny) {
    double t0 = 0.;
    ExactSol(Sol, x, y, t0, Nx, Ny);
}

/**
 * @brief Compute L2 norm of error array
 */
double TwoNorm(double** Err, int Nx, int Ny, double h) {
    auto it = array_map.find(Err);
    if (it == array_map.end()) {
        std::cerr << "ERROR [legacy]: TwoNorm called with unmapped array" << std::endl;
        return 0.0;
    }
    utils::Array2D* array = it->second;

    double e_2 = 0.0;
    for (int j = 1; j < Ny + 1; j++) {
        for (int i = 1; i < Nx + 1; i++) {
            e_2 += pow(Err[j][i], 2);
        }
    }
    return sqrt(e_2) * h;
}

/**
 * @brief Compute infinity norm of error array
 */
double InftyNorm(double** Err, int Nx, int Ny) {
    auto it = array_map.find(Err);
    if (it == array_map.end()) {
        std::cerr << "ERROR [legacy]: InftyNorm called with unmapped array" << std::endl;
        return 0.0;
    }
    utils::Array2D* array = it->second;

    double e_infty = 0.0;
    for (int j = 1; j < Ny + 1; j++) {
        for (int i = 1; i < Nx + 1; i++) {
            e_infty = std::max(e_infty, std::abs(Err[j][i]));
        }
    }
    return e_infty;
}

/**
 * @brief Compute error between numerical and reference solutions
 */
void Error(double** Sol, double** Ref, double** Err, int Nx, int Ny) {
    for (int j = 0; j < Ny + 2; j++) {
        for (int i = 0; i < Nx + 2; i++) {
            Err[j][i] = Sol[j][i] - Ref[j][i];
        }
    }
}

/**
 * @brief Copy contents of one 2D array to another
 */
void Copy(double** Sol, double** Sol0, int Nx, int Ny) {
    for (int j = 0; j < Ny + 2; j++) {
        for (int i = 0; i < Nx + 2; i++) {
            Sol0[j][i] = Sol[j][i];
        }
    }
}

/**
 * @brief Jacobi iteration solver for implicit Euler scheme (serial version)
 */
void Jacobi(double** x, double** x0, double** b, double& Residu, double Tol,
            int& iConv, int Nx, int Ny, double lambda) {
    double Coeff = 1.0 / (1.0 + 4.0 * lambda);
    double Res;
    Residu = 1.0;
    iConv = 0;

    while (Residu > Tol) {
        iConv = iConv + 1;
        Res = 0.0;
        for (int j = 1; j < Ny + 1; j++) {
            for (int i = 1; i < Nx + 1; i++) {
                // Jacobi update formula for 2D heat equation
                x[j][i] = Coeff * (lambda * (x0[j + 1][i] + x0[j][i + 1] +
                                              x0[j - 1][i] + x0[j][i - 1]) + b[j][i]);
                Res = std::max(std::abs(x[j][i] - x0[j][i]), Res);
            }
        }
        Residu = Res;
        Copy(x, x0, Nx, Ny);
    }
}

/**
 * @brief Export solution to file
 */
void Export(double** Sol, int I[], int J[], double x[], double y[],
            int Nx, int Ny, std::string SolFile) {
    std::ofstream outfile;
    outfile.open(SolFile);
    for (int j = 1; j < Ny + 1; j++) {
        for (int i = 1; i < Nx + 1; i++) {
            outfile << I[i] << " " << J[j] << " " << x[i] << " " << y[j] << " "
                    << std::setprecision(15) << Sol[j][i] << std::endl;
        }
    }
    outfile.close();
}

/**
 * @brief Cleanup function to free all mapped arrays
 * This should be called at the end of the program
 */
void legacy_cleanup() {
    for (auto& pair : array_map) {
        delete pair.second;
    }
    array_map.clear();
}
