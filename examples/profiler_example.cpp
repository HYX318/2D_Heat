/**
 * @file profiler_example.cpp
 * @brief Example demonstrating how to use the MPI Profiler class
 *
 * This example shows how to integrate the Profiler class into an MPI
 * application to measure communication and computation performance.
 */

#include "mpi/profiler.hpp"
#include "mpi/mpi_context.hpp"
#include <iostream>
#include <vector>
#include <cmath>

int main(int argc, char** argv) {
    // Initialize MPI
    MPIContext mpi(argc, argv);

    int rank = mpi.rank();
    int size = mpi.size();

    // Create profiler with PMPI enabled
    Profiler profiler(true, rank, size);

    if (rank == 0) {
        std::cout << "MPI Profiler Example\n";
        std::cout << "Running with " << size << " processes\n";
    }

    // Example 1: Simple ping-pong communication
    profiler.start_communication("ping_pong");

    std::vector<double> data(1000);
    std::fill(data.begin(), data.end(), rank);

    std::vector<double> recv_data(1000);
    int partner = (rank + 1) % size;

    // Record send
    profiler.record_send(partner, data.size() * sizeof(double), "ping_pong");

    MPI_Request send_req, recv_req;
    MPI_Isend(data.data(), data.size(), MPI_DOUBLE, partner, 0,
              MPI_COMM_WORLD, &send_req);
    MPI_Irecv_(recv_data.data(), recv_data.size(), MPI_DOUBLE,
               (rank + size - 1) % size, 0, MPI_COMM_WORLD, &recv_req);

    MPI_Wait(&send_req, MPI_STATUS_IGNORE);
    MPI_Wait(&recv_req, MPI_STATUS_IGNORE);

    // Record receive
    profiler.record_recv((rank + size - 1) % size,
                         recv_data.size() * sizeof(double), "ping_pong");

    profiler.end_communication("ping_pong");

    // Example 2: Computation phase
    profiler.start_computation("compute_step");

    std::vector<double> result(10000);
    for (size_t i = 0; i < result.size(); ++i) {
        result[i] = std::sin(i * 0.001) * std::cos(i * 0.002);
    }

    profiler.end_computation("compute_step");

    // Example 3: Global reduction with profiling
    profiler.start_communication("global_sum");

    double local_sum = std::accumulate(result.begin(), result.end(), 0.0);
    double global_sum = 0.0;

    profiler.record_send(0, sizeof(double), "sum_contrib");

    MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);

    if (rank == 0) {
        profiler.record_recv(MPI_ANY_SOURCE, sizeof(double), "sum_result");
    }

    profiler.end_communication("global_sum");

    // Synchronize and analyze load balance
    profiler.synchronize_and_analyze();

    // Print report on all processes
    profiler.print_report();

    // Also write to file
    std::string filename = "profiler_report_rank_" + std::to_string(rank) + ".txt";
    profiler.write_report(filename);

    // Write CSV and JSON for analysis
    std::string csv_file = "profiler_data_rank_" + std::to_string(rank) + ".csv";
    std::string json_file = "profiler_data_rank_" + std::to_string(rank) + ".json";

    profiler.write_csv(csv_file);
    profiler.write_json(json_file);

    if (rank == 0) {
        std::cout << "\nProfiler reports: " << filename << ", " << csv_file
                  << ", " << json_file << "\n";
    }

    return 0;
}
