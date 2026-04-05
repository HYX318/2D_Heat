/**
 * @file domain_decomposition_example.cpp
 * @brief Example demonstrating the use of DomainDecomposition class
 *
 * This example shows:
 * - Automatic domain decomposition
 * - Manual Cartesian decomposition
 * - Coordinate transformation
 * - Boundary detection
 * - Load balance analysis
 */

#include "mpi/mpi_context.hpp"
#include "mpi/cartesian_topology.hpp"
#include "mesh/domain_decomposition.hpp"
#include "utils/logger.hpp"
#include <iostream>
#include <iomanip>

int main(int argc, char** argv) {
    try {
        // Initialize MPI
        MPIContext mpi(argc, argv);
        Logger logger("DomainDecompositionExample");

        if (mpi.is_root()) {
            std::cout << "\n=== Domain Decomposition Example ===\n" << std::endl;
            std::cout << "Number of processes: " << mpi.size() << "\n" << std::endl;
        }

        // Create Cartesian topology with automatic dimensions
        CartesianTopology topology(MPI_COMM_WORLD);

        // ========================================
        // Example 1: Automatic decomposition
        // ========================================
        if (mpi.is_root()) {
            std::cout << "--- Example 1: Automatic Decomposition ---" << std::endl;
        }

        {
            size_t global_nx = 1000;
            size_t global_ny = 800;

            DomainDecomposition decomp(global_nx, global_ny, topology);

            const auto& my_domain = decomp.my_subdomain();

            std::cout << "Rank " << mpi.rank() << ": "
                      << "Domain [" << my_domain.start_i << ":"
                      << my_domain.start_i + my_domain.nx << "] x ["
                      << my_domain.start_j << ":"
                      << my_domain.start_j + my_domain.ny << "]"
                      << " (local size: " << my_domain.nx << "x" << my_domain.ny << ")"
                      << std::endl;

            // Coordinate transformation
            if (my_domain.nx > 0 && my_domain.ny > 0) {
                size_t global_i = my_domain.start_i + 5;
                size_t global_j = my_domain.start_j + 5;

                if (decomp.is_in_my_domain(global_i, global_j)) {
                    size_t local_i = my_domain.global_to_local_i(global_i);
                    size_t local_j = my_domain.global_to_local_j(global_j);
                    std::cout << "  Global (" << global_i << ", " << global_j << ") -> "
                              << "Local (" << local_i << ", " << local_j << ")" << std::endl;
                }
            }

            // Boundary detection
            std::cout << "  Boundaries: ";
            if (decomp.is_on_boundary(Direction::X)) {
                std::cout << "X ";
            }
            if (decomp.is_on_boundary(Direction::Y)) {
                std::cout << "Y ";
            }
            std::cout << std::endl;

            MPI_Barrier(MPI_COMM_WORLD);

            // Load balance analysis (only on root)
            if (mpi.is_root()) {
                auto stats = decomp.compute_load_balance();
                std::cout << "\nLoad Balance Statistics:" << std::endl;
                std::cout << "  Min load: " << stats.min_load << std::endl;
                std::cout << "  Max load: " << stats.max_load << std::endl;
                std::cout << "  Avg load: " << stats.avg_load << std::endl;
                std::cout << "  Imbalance ratio: " << std::fixed << std::setprecision(4)
                          << stats.imbalance_ratio << std::endl;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        // ========================================
        // Example 2: Manual Cartesian decomposition
        // ========================================
        if (mpi.is_root()) {
            std::cout << "\n--- Example 2: Manual Cartesian Decomposition ---" << std::endl;
        }

        if (mpi.size() >= 4) {
            size_t global_nx = 100;
            size_t global_ny = 100;
            size_t dim_x = 2;
            size_t dim_y = 2;

            DomainDecomposition decomp(global_nx, global_ny, dim_x, dim_y, topology);

            const auto& my_domain = decomp.my_subdomain();

            std::cout << "Rank " << mpi.rank() << ": "
                      << "Domain [" << my_domain.start_i << ":"
                      << my_domain.start_i + my_domain.nx << "] x ["
                      << my_domain.start_j << ":"
                      << my_domain.start_j + my_domain.ny << "]"
                      << std::endl;

            // Find owner for specific coordinates
            if (mpi.is_root()) {
                int owner = decomp.find_owner_rank(25, 25);
                std::cout << "Cell (25, 25) is owned by rank " << owner << std::endl;
            }
        } else if (mpi.is_root()) {
            std::cout << "Skipped: Need at least 4 processes" << std::endl;
        }

        MPI_Barrier(MPI_COMM_WORLD);

        // ========================================
        // Example 3: Manual non-uniform decomposition
        // ========================================
        if (mpi.is_root()) {
            std::cout << "\n--- Example 3: Manual Non-Uniform Decomposition ---" << std::endl;
        }

        if (mpi.size() == 2) {
            size_t global_nx = 100;
            size_t global_ny = 100;

            // Create non-uniform decomposition
            std::vector<Subdomain> subdomains(2);
            subdomains[0] = {0, 0, 75, 100, global_nx, global_ny};  // Rank 0: 75x100
            subdomains[1] = {75, 0, 25, 100, global_nx, global_ny}; // Rank 1: 25x100

            DomainDecomposition decomp(subdomains, topology);

            const auto& my_domain = decomp.my_subdomain();

            std::cout << "Rank " << mpi.rank() << ": "
                      << "Domain [" << my_domain.start_i << ":"
                      << my_domain.start_i + my_domain.nx << "] x ["
                      << my_domain.start_j << ":"
                      << my_domain.start_j + my_domain.ny << "]"
                      << std::endl;

            MPI_Barrier(MPI_COMM_WORLD);

            // Load balance analysis
            if (mpi.is_root()) {
                auto stats = decomp.compute_load_balance();
                std::cout << "Load Balance Statistics:" << std::endl;
                std::cout << "  Min load: " << stats.min_load << std::endl;
                std::cout << "  Max load: " << stats.max_load << std::endl;
                std::cout << "  Avg load: " << stats.avg_load << std::endl;
                std::cout << "  Imbalance ratio: " << std::fixed << std::setprecision(4)
                          << stats.imbalance_ratio << std::endl;
            }
        } else if (mpi.is_root()) {
            std::cout << "Skipped: Need exactly 2 processes" << std::endl;
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (mpi.is_root()) {
            std::cout << "\n=== Example Complete ===" << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
