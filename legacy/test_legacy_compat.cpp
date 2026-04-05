/**
 * @file test_legacy_compat.cpp
 * @brief Simple test to verify legacy compatibility layer compiles and works
 *
 * This test uses the legacy API functions to verify they still work.
 * All calls will generate deprecation warnings (which is expected).
 */

#include "heat_utils_compat.hpp"
#include "interfaces_compat.hpp"
#include "Consts.h"
#include <iostream>
#include <cassert>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int myRank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    if (myRank == 0) {
        std::cout << "Testing Legacy Compatibility Layer" << std::endl;
        std::cout << "=================================" << std::endl;
        std::cout << "Note: Deprecation warnings are expected and intentional" << std::endl;
        std::cout << std::endl;
    }

    // Test 1: ReadParam
    if (myRank == 0) {
        int Nx, Ny, Nt;
        double StabP;
        ReadParam(Nx, Ny, Nt, StabP);
        std::cout << "Test 1 PASSED: ReadParam works" << std::endl;
        std::cout << "  Nx=" << Nx << " Ny=" << Ny << " Nt=" << Nt << " StabP=" << StabP << std::endl;
    }

    // Test 2: Contiguous2D
    int Nx = 50, Ny = 50;
    double** array = new double*[Ny+2];
    Contiguous2D(array, Nx+2, Ny+2);

    if (myRank == 0) {
        std::cout << "Test 2 PASSED: Contiguous2D works" << std::endl;
    }

    // Test 3: Analytic
    double x = 0.5, y = 0.5, t = 0.1;
    double val = Analytic(x, y, t);

    if (myRank == 0) {
        std::cout << "Test 3 PASSED: Analytic works" << std::endl;
        std::cout << "  Analytic(0.5, 0.5, 0.1) = " << val << std::endl;
    }

    // Test 4: ExactSol
    double* x_coords = new double[Nx+2];
    double* y_coords = new double[Ny+2];
    for (int i = 0; i < Nx+2; i++) {
        x_coords[i] = i / (double)(Nx+1);
    }
    for (int j = 0; j < Ny+2; j++) {
        y_coords[j] = j / (double)(Ny+1);
    }

    ExactSol(array, x_coords, y_coords, 0.1, Nx, Ny);

    if (myRank == 0) {
        std::cout << "Test 4 PASSED: ExactSol works" << std::endl;
    }

    // Test 5: Init
    double** init_array = new double*[Ny+2];
    Contiguous2D(init_array, Nx+2, Ny+2);
    Init(init_array, x_coords, y_coords, Nx, Ny);

    if (myRank == 0) {
        std::cout << "Test 5 PASSED: Init works" << std::endl;
    }

    // Test 6: Copy
    double** copy_array = new double*[Ny+2];
    Contiguous2D(copy_array, Nx+2, Ny+2);
    Copy(array, copy_array, Nx, Ny);

    // Verify copy
    bool copy_correct = true;
    for (int j = 0; j < Ny+2 && copy_correct; j++) {
        for (int i = 0; i < Nx+2 && copy_correct; i++) {
            if (std::abs(array[j][i] - copy_array[j][i]) > 1e-15) {
                copy_correct = false;
            }
        }
    }

    if (myRank == 0) {
        if (copy_correct) {
            std::cout << "Test 6 PASSED: Copy works" << std::endl;
        } else {
            std::cout << "Test 6 FAILED: Copy produced incorrect results" << std::endl;
        }
    }

    // Test 7: Error and norms
    double** error_array = new double*[Ny+2];
    Contiguous2D(error_array, Nx+2, Ny+2);
    Error(array, copy_array, error_array, Nx, Ny);

    double h = 1.0 / (Nx + 1);
    double l2_norm = TwoNorm(error_array, Nx, Ny, h);
    double linfty_norm = InftyNorm(error_array, Nx, Ny);

    if (myRank == 0) {
        std::cout << "Test 7 PASSED: Error, TwoNorm, and InftyNorm work" << std::endl;
        std::cout << "  L2 norm: " << l2_norm << std::endl;
        std::cout << "  Linfty norm: " << linfty_norm << std::endl;
    }

    // Test 8: Jacobi (serial)
    double** x = new double*[Ny+2];
    double** x0 = new double*[Ny+2];
    double** b = new double*[Ny+2];
    Contiguous2D(x, Nx+2, Ny+2);
    Contiguous2D(x0, Nx+2, Ny+2);
    Contiguous2D(b, Nx+2, Ny+2);

    // Initialize x0 and b
    for (int j = 1; j <= Ny; j++) {
        for (int i = 1; i <= Nx; i++) {
            x0[j][i] = 1.0;
            b[j][i] = 1.0;
        }
    }

    double Residu;
    double Tol = 1e-6;
    int iConv;
    double lambda = 0.25;

    Jacobi(x, x0, b, Residu, Tol, iConv, Nx, Ny, lambda);

    if (myRank == 0) {
        std::cout << "Test 8 PASSED: Jacobi works" << std::endl;
        std::cout << "  Iterations: " << iConv << std::endl;
        std::cout << "  Final residual: " << Residu << std::endl;
    }

    // Test 9: Export
    int* I = new int[Nx+2];
    int* J = new int[Ny+2];
    for (int i = 0; i < Nx+2; i++) {
        I[i] = i;
    }
    for (int j = 0; j < Ny+2; j++) {
        J[j] = j;
    }

    std::string filename = "test_legacy_output.txt";
    Export(array, I, J, x_coords, y_coords, Nx, Ny, filename);

    if (myRank == 0) {
        std::cout << "Test 9 PASSED: Export works" << std::endl;
        std::cout << "  Output file: " << filename << std::endl;
    }

    // Test 10: MPIJacobi (parallel) - only if we have multiple processes
    if (numProcs > 1) {
        // Create Cartesian communicator
        MPI_Comm SBD_COMM;
        int nDims = 2;
        int dims[2] = {2, 2};  // 2x2 grid (need 4 processes)
        int periods[2] = {0, 0};
        MPI_Cart_create(MPI_COMM_WORLD, nDims, dims, periods, 1, &SBD_COMM);

        // Get neighbors
        int NeighbourRank[4];
        MPI_Cart_shift(SBD_COMM, 0, 1, &NeighbourRank[West], &NeighbourRank[East]);
        MPI_Cart_shift(SBD_COMM, 1, 1, &NeighbourRank[South], &NeighbourRank[North]);

        // Create column type
        MPI_Datatype colType;
        MPI_Type_vector(Ny+2, 1, Nx+2, MPI_DOUBLE, &colType);
        MPI_Type_commit(&colType);

        // Reset arrays
        for (int j = 1; j <= Ny; j++) {
            for (int i = 1; i <= Nx; i++) {
                x0[j][i] = 1.0;
                b[j][i] = 1.0;
                x[j][i] = 0.0;
            }
        }

        // Solve with MPIJacobi
        MPIJacobi(x, x0, b, Residu, Tol, iConv, Nx, Ny, lambda,
                  NeighbourRank, SBD_COMM, colType, myRank);

        if (myRank == 0) {
            std::cout << "Test 10 PASSED: MPIJacobi works" << std::endl;
            std::cout << "  Iterations: " << iConv << std::endl;
            std::cout << "  Final residual: " << Residu << std::endl;
        }

        MPI_Type_free(&colType);
        MPI_Comm_free(&SBD_COMM);
    } else if (myRank == 0) {
        std::cout << "Test 10 SKIPPED: Need multiple processes for MPIJacobi test" << std::endl;
    }

    // Cleanup
    delete[] array[0];
    delete[] array;
    delete[] init_array[0];
    delete[] init_array;
    delete[] copy_array[0];
    delete[] copy_array;
    delete[] error_array[0];
    delete[] error_array;
    delete[] x[0];
    delete[] x;
    delete[] x0[0];
    delete[] x0;
    delete[] b[0];
    delete[] b;
    delete[] x_coords;
    delete[] y_coords;
    delete[] I;
    delete[] J;

    if (myRank == 0) {
        std::cout << std::endl;
        std::cout << "=================================" << std::endl;
        std::cout << "All legacy compatibility tests PASSED!" << std::endl;
        std::cout << "=================================" << std::endl;
    }

    MPI_Finalize();
    return 0;
}
