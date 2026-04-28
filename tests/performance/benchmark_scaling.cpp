/**
 * MPI Scaling Performance Benchmark
 *
 * This file implements scaling tests for MPI applications:
 * - Strong scaling: Fixed total workload, increasing processors
 * - Weak scaling: Fixed workload per processor, increasing processors
 */

#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <numeric>

using namespace std;

// Helper class for CSV output
class CSVWriter {
public:
    CSVWriter(int rank, const std::string& filename)
        : rank_(rank), filename_(filename) {}

    void write_header() {
        if (rank_ == 0) {
            ofstream file(filename_);
            if (file.is_open()) {
                file << "test_type,processes,problem_size,per_proc_size,time_ms,efficiency,gflops\n";
                file.close();
            }
        }
    }

    void write(const std::string& test_type, int processes, size_t problem_size,
               size_t per_proc_size, double time_ms, double efficiency, double gflops) {
        if (rank_ == 0) {
            ofstream file(filename_, ios::app);
            if (file.is_open()) {
                file << test_type << "," << processes << "," << problem_size << ","
                     << per_proc_size << "," << fixed << setprecision(4) << time_ms << ","
                     << fixed << setprecision(4) << efficiency << ","
                     << fixed << setprecision(4) << gflops << "\n";
                file.close();
            }
        }
    }

private:
    int rank_;
    std::string filename_;
};

// Simulate computation on a 2D grid
void simulate_computation(size_t grid_size, int iterations) {
    vector<double> grid(grid_size * grid_size, 0.0);
    vector<double> new_grid(grid_size * grid_size, 0.0);

    for (int iter = 0; iter < iterations; iter++) {
        for (size_t i = 1; i < grid_size - 1; i++) {
            for (size_t j = 1; j < grid_size - 1; j++) {
                size_t idx = i * grid_size + j;
                new_grid[idx] = 0.25 * (grid[idx - 1] + grid[idx + 1] +
                                        grid[idx - grid_size] + grid[idx + grid_size]);
            }
        }
        swap(grid, new_grid);
    }
}

// Simulate communication between neighboring processes
double simulate_communication(int rank, int size, int iterations, size_t message_size) {
    vector<double> send_buffer(message_size, static_cast<double>(rank));
    vector<double> recv_buffer(message_size, 0.0);

    double start_time = MPI_Wtime();

    for (int i = 0; i < iterations; i++) {
        // Communicate with neighbors
        int left = (rank > 0) ? rank - 1 : MPI_PROC_NULL;
        int right = (rank < size - 1) ? rank + 1 : MPI_PROC_NULL;

        MPI_Request send_req[2], recv_req[2];
        int req_count = 0;

        // Non-blocking sends
        if (right != MPI_PROC_NULL) {
            MPI_Isend(send_buffer.data(), message_size, MPI_DOUBLE, right, 0,
                     MPI_COMM_WORLD, &send_req[req_count]);
            req_count++;
        }
        if (left != MPI_PROC_NULL) {
            MPI_Isend(send_buffer.data(), message_size, MPI_DOUBLE, left, 0,
                     MPI_COMM_WORLD, &send_req[req_count]);
            req_count++;
        }

        req_count = 0;

        // Non-blocking receives
        if (left != MPI_PROC_NULL) {
            MPI_Irecv(recv_buffer.data(), message_size, MPI_DOUBLE, left, 0,
                     MPI_COMM_WORLD, &recv_req[req_count]);
            req_count++;
        }
        if (right != MPI_PROC_NULL) {
            MPI_Irecv(recv_buffer.data(), message_size, MPI_DOUBLE, right, 0,
                     MPI_COMM_WORLD, &recv_req[req_count]);
            req_count++;
        }

        // Wait for all communications to complete
        MPI_Waitall(2, recv_req, MPI_STATUSES_IGNORE);
        MPI_Waitall(2, send_req, MPI_STATUSES_IGNORE);
    }

    return MPI_Wtime() - start_time;
}

// Strong scaling test: Fixed total workload, increasing processors
void test_strong_scaling(int rank, int size, size_t total_grid_size,
                         int compute_iterations, int comm_iterations,
                         CSVWriter& csv) {
    if (rank == 0) {
        cout << "\n=== Strong Scaling Test ===\n";
        cout << "Total grid size: " << total_grid_size << " x " << total_grid_size << "\n";
        cout << "Compute iterations: " << compute_iterations << "\n";
        cout << "Comm iterations: " << comm_iterations << "\n";
        cout << "Processes: " << size << "\n\n";
    }

    // Calculate per-process grid size
    size_t per_proc_grid = max(static_cast<size_t>(1), total_grid_size / static_cast<size_t>(size));

    // Warm-up
    simulate_computation(per_proc_grid, 5);
    simulate_communication(rank, size, 5, per_proc_grid);

    MPI_Barrier(MPI_COMM_WORLD);

    // Benchmark computation
    double compute_start = MPI_Wtime();
    simulate_computation(per_proc_grid, compute_iterations);
    double compute_time = (MPI_Wtime() - compute_start) * 1000.0; // Convert to ms

    // Benchmark communication
    double comm_time = simulate_communication(rank, size, comm_iterations, per_proc_grid) * 1000.0;

    // Calculate total time
    double total_time = compute_time + comm_time;

    // Calculate efficiency (theoretical speedup = size)
    // Efficiency = (T1 / (Tp * p)) * 100
    // For strong scaling, we compare against a hypothetical single processor time
    double baseline_time = compute_time * size + comm_time; // Approximate baseline
    double speedup = baseline_time / total_time;
    double efficiency = (speedup / size) * 100.0;

    // Calculate GFLOPS (rough estimate for stencil computation)
    // 5 floating point operations per cell per iteration
    double flops = static_cast<double>(per_proc_grid) * static_cast<double>(per_proc_grid) *
                   compute_iterations * 5.0 * size;
    double gflops = flops / (total_time * 1e-6);

    // Reduce to get maximum time across all processes
    double max_compute_time, max_comm_time, max_total_time;
    MPI_Reduce(&compute_time, &max_compute_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&comm_time, &max_comm_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_time, &max_total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "Results:\n";
        cout << "  Compute time:  " << fixed << setprecision(2) << max_compute_time << " ms\n";
        cout << "  Communication time: " << fixed << setprecision(2) << max_comm_time << " ms\n";
        cout << "  Total time:    " << fixed << setprecision(2) << max_total_time << " ms\n";
        cout << "  Efficiency:    " << fixed << setprecision(2) << efficiency << "%\n";
        cout << "  Performance:  " << fixed << setprecision(2) << gflops << " GFLOPS\n";

        csv.write("strong_scaling", size, total_grid_size * total_grid_size,
                  per_proc_grid * per_proc_grid, max_total_time, efficiency, gflops);
    }
}

// Weak scaling test: Fixed workload per processor, increasing processors
void test_weak_scaling(int rank, int size, size_t per_proc_grid_size,
                       int compute_iterations, int comm_iterations,
                       CSVWriter& csv) {
    if (rank == 0) {
        cout << "\n=== Weak Scaling Test ===\n";
        cout << "Per-process grid size: " << per_proc_grid_size << " x " << per_proc_grid_size << "\n";
        cout << "Compute iterations: " << compute_iterations << "\n";
        cout << "Comm iterations: " << comm_iterations << "\n";
        cout << "Processes: " << size << "\n\n";
    }

    // Total problem size scales with number of processes
    size_t total_grid_size = per_proc_grid_size * static_cast<size_t>(size);

    // Warm-up
    simulate_computation(per_proc_grid_size, 5);
    simulate_communication(rank, size, 5, per_proc_grid_size);

    MPI_Barrier(MPI_COMM_WORLD);

    // Benchmark computation
    double compute_start = MPI_Wtime();
    simulate_computation(per_proc_grid_size, compute_iterations);
    double compute_time = (MPI_Wtime() - compute_start) * 1000.0; // Convert to ms

    // Benchmark communication
    double comm_time = simulate_communication(rank, size, comm_iterations, per_proc_grid_size) * 1000.0;

    // Calculate total time
    double total_time = compute_time + comm_time;

    // For weak scaling, ideal efficiency is 100% (constant time per process)
    // We calculate efficiency based on how close we are to ideal
    // Ideal: time per process remains constant
    // Efficiency = T1 / Tp * 100 (where T1 is time for single process)
    // Since we can't easily measure single process time, we use theoretical efficiency
    // For weak scaling, efficiency should be ~100% if communication scales well
    double comm_ratio = comm_time / total_time;
    double efficiency = 100.0 * (1.0 - comm_ratio); // Simplified efficiency metric

    // Calculate total GFLOPS
    double flops = static_cast<double>(per_proc_grid_size) * static_cast<double>(per_proc_grid_size) *
                   compute_iterations * 5.0 * size;
    double gflops = flops / (total_time * 1e-6);

    // Reduce to get maximum time across all processes
    double max_compute_time, max_comm_time, max_total_time;
    MPI_Reduce(&compute_time, &max_compute_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&comm_time, &max_comm_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_time, &max_total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "Results:\n";
        cout << "  Total grid size: " << total_grid_size << "\n";
        cout << "  Compute time:  " << fixed << setprecision(2) << max_compute_time << " ms\n";
        cout << "  Communication time: " << fixed << setprecision(2) << max_comm_time << " ms\n";
        cout << "  Total time:    " << fixed << setprecision(2) << max_total_time << " ms\n";
        cout << "  Efficiency:    " << fixed << setprecision(2) << efficiency << "%\n";
        cout << "  Performance:  " << fixed << setprecision(2) << gflops << " GFLOPS\n";

        csv.write("weak_scaling", size, total_grid_size * total_grid_size,
                  per_proc_grid_size * per_proc_grid_size, max_total_time, efficiency, gflops);
    }
}

// Run scaling test across different process counts
void run_scaling_suite(int rank, int size, int max_procs,
                       size_t strong_grid_size, size_t weak_per_proc_size,
                       int compute_iterations, int comm_iterations,
                       CSVWriter& csv) {
    // Only run at the maximum number of processes
    // The test suite should be run multiple times with different process counts

    if (size == max_procs) {
        if (rank == 0) {
            cout << "\n========================================\n";
            cout << "   Scaling Test Suite\n";
            cout << "========================================\n";
            cout << "Maximum processes: " << max_procs << "\n";
            cout << "Note: Run this test with different process counts\n";
            cout << "      to generate scaling curves:\n";
            cout << "      mpirun -np 1  ./benchmark_scaling ...\n";
            cout << "      mpirun -np 2  ./benchmark_scaling ...\n";
            cout << "      mpirun -np 4  ./benchmark_scaling ...\n";
            cout << "      mpirun -np 8  ./benchmark_scaling ...\n";
            cout << "      mpirun -np 16 ./benchmark_scaling ...\n";
        }
    }

    // Run strong scaling
    test_strong_scaling(rank, size, strong_grid_size, compute_iterations, comm_iterations, csv);

    // Run weak scaling
    test_weak_scaling(rank, size, weak_per_proc_size, compute_iterations, comm_iterations, csv);
}

// Print usage information
void print_usage(const char* prog_name) {
    cout << "Usage: " << prog_name << " [options]\n\n";
    cout << "Options:\n";
    cout << "  --type <type>            Test type: strong, weak, or both (default: both)\n";
    cout << "  --max <n>                Maximum number of processes (default: current)\n";
    cout << "  --strong-size <n>        Total grid size for strong scaling (default: 1024)\n";
    cout << "  --weak-per-proc <n>       Per-process grid size for weak scaling (default: 256)\n";
    cout << "  --compute-iter <n>       Number of compute iterations (default: 100)\n";
    cout << "  --comm-iter <n>          Number of communication iterations (default: 1000)\n";
    cout << "  --output <file>          CSV output filename (default: scaling_results.csv)\n";
    cout << "  --help                   Show this help message\n";
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Default parameters
    string test_type = "both";
    int max_procs = size;
    size_t strong_grid_size = 1024;
    size_t weak_per_proc_size = 256;
    int compute_iterations = 100;
    int comm_iterations = 1000;
    string csv_filename = "scaling_results.csv";

    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];

        if (arg == "--help") {
            if (rank == 0) {
                print_usage(argv[0]);
            }
            MPI_Finalize();
            return 0;
        } else if (arg == "--type") {
            if (i + 1 < argc) {
                test_type = argv[++i];
            }
        } else if (arg == "--max") {
            if (i + 1 < argc) {
                max_procs = atoi(argv[++i]);
            }
        } else if (arg == "--strong-size") {
            if (i + 1 < argc) {
                strong_grid_size = atoi(argv[++i]);
            }
        } else if (arg == "--weak-per-proc") {
            if (i + 1 < argc) {
                weak_per_proc_size = atoi(argv[++i]);
            }
        } else if (arg == "--compute-iter") {
            if (i + 1 < argc) {
                compute_iterations = atoi(argv[++i]);
            }
        } else if (arg == "--comm-iter") {
            if (i + 1 < argc) {
                comm_iterations = atoi(argv[++i]);
            }
        } else if (arg == "--output") {
            if (i + 1 < argc) {
                csv_filename = argv[++i];
            }
        }
    }

    CSVWriter csv(rank, csv_filename);

    if (rank == 0) {
        cout << "\n========================================\n";
        cout << "   MPI Scaling Performance Benchmark\n";
        cout << "========================================\n";
        cout << "Current processes: " << size << "\n";
        cout << "CSV output: " << csv_filename << "\n";
    }

    csv.write_header();

    // Run the scaling suite
    run_scaling_suite(rank, size, max_procs, strong_grid_size, weak_per_proc_size,
                       compute_iterations, comm_iterations, csv);

    if (rank == 0) {
        cout << "\n========================================\n";
        cout << "   Scaling Test Complete\n";
        cout << "========================================\n";
        cout << "Results saved to: " << csv_filename << "\n";
    }

    MPI_Finalize();
    return 0;
}
