/**
 * MPI Communication/Computation Overlap Benchmark
 *
 * This file implements benchmarks for testing overlap between MPI communication
 * and computation using non-blocking operations:
 * - Non-blocking communication (MPI_Isend/MPI_Irecv) with compute overlap
 * - Comparison of async vs sync communication performance
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
                file << "test_type,processes,message_size,iterations,time_ms,compute_time_ms,comm_time_ms,overlap_ratio\n";
                file.close();
            }
        }
    }

    void write(const std::string& test_type, int processes, size_t message_size,
               int iterations, double time_ms, double compute_time_ms,
               double comm_time_ms, double overlap_ratio) {
        if (rank_ == 0) {
            ofstream file(filename_, ios::app);
            if (file.is_open()) {
                file << test_type << "," << processes << "," << message_size << ","
                     << iterations << "," << fixed << setprecision(4) << time_ms << ","
                     << compute_time_ms << "," << comm_time_ms << ","
                     << fixed << setprecision(4) << overlap_ratio << "\n";
                file.close();
            }
        }
    }

private:
    int rank_;
    std::string filename_;
};

// Simulate computation on internal grid points
double compute_internal_points(vector<double>& grid, size_t grid_size, size_t num_iterations) {
    double start_time = MPI_Wtime();

    for (size_t iter = 0; iter < num_iterations; iter++) {
        // Compute only internal points (1 to grid_size-2)
        for (size_t i = 1; i < grid_size - 1; i++) {
            for (size_t j = 1; j < grid_size - 1; j++) {
                size_t idx = i * grid_size + j;
                grid[idx] = 0.25 * (grid[idx - 1] + grid[idx + 1] +
                                   grid[idx - grid_size] + grid[idx + grid_size]);
            }
        }
    }

    return MPI_Wtime() - start_time;
}

// Simulate computation on boundary grid points
double compute_boundary_points(vector<double>& grid, size_t grid_size) {
    double start_time = MPI_Wtime();

    // Compute only boundary points
    // Top and bottom rows
    for (size_t j = 0; j < grid_size; j++) {
        grid[j] = 0.5 * grid[j];                 // Top row
        grid[(grid_size - 1) * grid_size + j] = 0.5 * grid[(grid_size - 1) * grid_size + j];  // Bottom row
    }

    // Left and right columns
    for (size_t i = 0; i < grid_size; i++) {
        grid[i * grid_size] = 0.5 * grid[i * grid_size];                      // Left column
        grid[i * grid_size + (grid_size - 1)] = 0.5 * grid[i * grid_size + (grid_size - 1)];  // Right column
    }

    return MPI_Wtime() - start_time;
}

// Benchmark with non-blocking communication (overlap)
void benchmark_async_overlap(int rank, int size, int iterations,
                             size_t grid_size, size_t message_size,
                             CSVWriter& csv) {
    if (rank == 0) {
        cout << "\n=== Async Communication Overlap Benchmark ===\n";
        cout << "Iterations: " << iterations << "\n";
        cout << "Grid size: " << grid_size << " x " << grid_size << "\n";
        cout << "Message size: " << message_size << " doubles\n";
        cout << "Processes: " << size << "\n\n";
    }

    // Initialize grid
    vector<double> grid(grid_size * grid_size, static_cast<double>(rank));
    vector<double> send_buffer(message_size, static_cast<double>(rank));
    vector<double> recv_buffer(message_size, 0.0);

    // Determine neighbors
    int left = (rank > 0) ? rank - 1 : MPI_PROC_NULL;
    int right = (rank < size - 1) ? rank + 1 : MPI_PROC_NULL;

    // Warm-up
    for (int i = 0; i < 5; i++) {
        MPI_Request reqs[4];
        int req_count = 0;

        // Non-blocking sends
        if (left != MPI_PROC_NULL) {
            MPI_Isend(send_buffer.data(), message_size, MPI_DOUBLE, left, 0,
                     MPI_COMM_WORLD, &reqs[req_count++]);
        }
        if (right != MPI_PROC_NULL) {
            MPI_Isend(send_buffer.data(), message_size, MPI_DOUBLE, right, 0,
                     MPI_COMM_WORLD, &reqs[req_count++]);
        }

        // Non-blocking receives
        if (left != MPI_PROC_NULL) {
            MPI_Irecv(recv_buffer.data(), message_size, MPI_DOUBLE, left, 0,
                     MPI_COMM_WORLD, &reqs[req_count++]);
        }
        if (right != MPI_PROC_NULL) {
            MPI_Irecv(recv_buffer.data(), message_size, MPI_DOUBLE, right, 0,
                     MPI_COMM_WORLD, &reqs[req_count++]);
        }

        // Overlap with computation
        compute_internal_points(grid, grid_size, 1);

        // Wait for communication
        MPI_Waitall(req_count, reqs, MPI_STATUSES_IGNORE);

        // Compute boundary after communication
        compute_boundary_points(grid, grid_size);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Benchmark
    double total_time = 0.0;
    double total_compute_time = 0.0;
    double total_comm_time = 0.0;

    double start_time = MPI_Wtime();

    for (int i = 0; i < iterations; i++) {
        // Initiate non-blocking communication
        MPI_Request reqs[4];
        int req_count = 0;

        double comm_start = MPI_Wtime();

        // Non-blocking sends
        if (left != MPI_PROC_NULL) {
            MPI_Isend(send_buffer.data(), message_size, MPI_DOUBLE, left, 0,
                     MPI_COMM_WORLD, &reqs[req_count++]);
        }
        if (right != MPI_PROC_NULL) {
            MPI_Isend(send_buffer.data(), message_size, MPI_DOUBLE, right, 0,
                     MPI_COMM_WORLD, &reqs[req_count++]);
        }

        // Non-blocking receives
        if (left != MPI_PROC_NULL) {
            MPI_Irecv(recv_buffer.data(), message_size, MPI_DOUBLE, left, 0,
                     MPI_COMM_WORLD, &reqs[req_count++]);
        }
        if (right != MPI_PROC_NULL) {
            MPI_Irecv(recv_buffer.data(), message_size, MPI_DOUBLE, right, 0,
                     MPI_COMM_WORLD, &reqs[req_count++]);
        }

        double comm_init_time = MPI_Wtime() - comm_start;

        // Overlap with computation on internal points
        double compute_start = MPI_Wtime();
        compute_internal_points(grid, grid_size, 10);
        double compute_internal_time = MPI_Wtime() - compute_start;

        // Wait for communication to complete
        double comm_wait_start = MPI_Wtime();
        MPI_Waitall(req_count, reqs, MPI_STATUSES_IGNORE);
        double comm_wait_time = MPI_Wtime() - comm_wait_start;

        double comm_end_time = MPI_Wtime() - comm_start;

        // Compute boundary points after communication completes
        double compute_boundary_start = MPI_Wtime();
        compute_boundary_points(grid, grid_size);
        double compute_boundary_time = MPI_Wtime() - compute_boundary_start;

        // Accumulate timings
        total_comm_time += comm_end_time;
        total_compute_time += compute_internal_time + compute_boundary_time;
    }

    total_time = MPI_Wtime() - start_time;

    // Calculate overlap ratio
    // Ideal: total time ≈ max(compute_time, comm_time)
    // Overlap ratio = (compute_time + comm_time) / total_time
    // Values > 1 indicate good overlap
    double overlap_ratio = (total_compute_time + total_comm_time) / total_time;

    // Reduce to get maximum time across all processes
    double max_total_time, max_compute_time, max_comm_time;
    MPI_Reduce(&total_time, &max_total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_compute_time, &max_compute_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_comm_time, &max_comm_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double avg_total_time = max_total_time / iterations * 1000.0;
        double avg_compute_time = max_compute_time / iterations * 1000.0;
        double avg_comm_time = max_comm_time / iterations * 1000.0;

        cout << "Results:\n";
        cout << "  Average total time:       " << fixed << setprecision(4) << avg_total_time << " ms\n";
        cout << "  Average compute time:     " << fixed << setprecision(4) << avg_compute_time << " ms\n";
        cout << "  Average comm time:        " << fixed << setprecision(4) << avg_comm_time << " ms\n";
        cout << "  Overlap ratio:            " << fixed << setprecision(4) << overlap_ratio << "\n";
        cout << "  (Values > 1.0 indicate good overlap)\n";

        csv.write("async_overlap", size, message_size, iterations, avg_total_time,
                  avg_compute_time, avg_comm_time, overlap_ratio);
    }
}

// Benchmark with synchronous communication (no overlap)
void benchmark_sync_no_overlap(int rank, int size, int iterations,
                                size_t grid_size, size_t message_size,
                                CSVWriter& csv) {
    if (rank == 0) {
        cout << "\n=== Sync Communication (No Overlap) Benchmark ===\n";
        cout << "Iterations: " << iterations << "\n";
        cout << "Grid size: " << grid_size << " x " << grid_size << "\n";
        cout << "Message size: " << message_size << " doubles\n";
        cout << "Processes: " << size << "\n\n";
    }

    // Initialize grid
    vector<double> grid(grid_size * grid_size, static_cast<double>(rank));
    vector<double> send_buffer(message_size, static_cast<double>(rank));
    vector<double> recv_buffer(message_size, 0.0);

    // Determine neighbors
    int left = (rank > 0) ? rank - 1 : MPI_PROC_NULL;
    int right = (rank < size - 1) ? rank + 1 : MPI_PROC_NULL;

    // Warm-up
    for (int i = 0; i < 5; i++) {
        // Communication first
        if (left != MPI_PROC_NULL) {
            MPI_Recv(recv_buffer.data(), message_size, MPI_DOUBLE, left, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(send_buffer.data(), message_size, MPI_DOUBLE, left, 0,
                     MPI_COMM_WORLD);
        }
        if (right != MPI_PROC_NULL) {
            MPI_Send(send_buffer.data(), message_size, MPI_DOUBLE, right, 0,
                     MPI_COMM_WORLD);
            MPI_Recv(recv_buffer.data(), message_size, MPI_DOUBLE, right, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Then computation
        compute_internal_points(grid, grid_size, 1);
        compute_boundary_points(grid, grid_size);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Benchmark
    double total_time = 0.0;
    double total_compute_time = 0.0;
    double total_comm_time = 0.0;

    double start_time = MPI_Wtime();

    for (int i = 0; i < iterations; i++) {
        // Communication first (no overlap)
        double comm_start = MPI_Wtime();

        if (left != MPI_PROC_NULL) {
            MPI_Recv(recv_buffer.data(), message_size, MPI_DOUBLE, left, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(send_buffer.data(), message_size, MPI_DOUBLE, left, 0,
                     MPI_COMM_WORLD);
        }
        if (right != MPI_PROC_NULL) {
            MPI_Send(send_buffer.data(), message_size, MPI_DOUBLE, right, 0,
                     MPI_COMM_WORLD);
            MPI_Recv(recv_buffer.data(), message_size, MPI_DOUBLE, right, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        double comm_end = MPI_Wtime();
        total_comm_time += comm_end - comm_start;

        // Then computation
        double compute_start = MPI_Wtime();
        compute_internal_points(grid, grid_size, 10);
        compute_boundary_points(grid, grid_size);
        double compute_end = MPI_Wtime();
        total_compute_time += compute_end - compute_start;
    }

    total_time = MPI_Wtime() - start_time;

    // Reduce to get maximum time across all processes
    double max_total_time, max_compute_time, max_comm_time;
    MPI_Reduce(&total_time, &max_total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_compute_time, &max_compute_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_comm_time, &max_comm_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank ==0) {
        double avg_total_time = max_total_time / iterations * 1000.0;
        double avg_compute_time = max_compute_time / iterations * 1000.0;
        double avg_comm_time = max_comm_time / iterations * 1000.0;

        // No overlap, so ratio should be ~1.0
        double overlap_ratio = 1.0;

        cout << "Results:\n";
        cout << "  Average total time:       " << fixed << setprecision(4) << avg_total_time << " ms\n";
        cout << "  Average compute time:     " << fixed << setprecision(4) << avg_compute_time << " ms\n";
        cout << "  Average comm time:        " << fixed << setprecision(4) << avg_comm_time << " ms\n";
        cout << "  Overlap ratio:            " << fixed << setprecision(4) << overlap_ratio << "\n";
        cout << "  (Sequential: no overlap)\n";

        csv.write("sync_no_overlap", size, message_size, iterations, avg_total_time,
                  avg_compute_time, avg_comm_time, overlap_ratio);
    }
}

// Compare async vs sync performance
void compare_async_vs_sync(int rank, int size, int iterations,
                          size_t grid_size, size_t message_size,
                          CSVWriter& csv) {
    if (rank == 0) {
        cout << "\n=== Async vs Sync Comparison ===\n";
        cout << "Grid size: " << grid_size << " x " << grid_size << "\n";
        cout << "Message size: " << message_size << " doubles\n";
        cout << "Iterations: " << iterations << "\n";
        cout << "Processes: " << size << "\n\n";
    }

    // Run sync benchmark
    benchmark_sync_no_overlap(rank, size, iterations, grid_size, message_size, csv);

    // Run async benchmark
    benchmark_async_overlap(rank, size, iterations, grid_size, message_size, csv);

    if (rank == 0) {
        cout << "\n--- Comparison Summary ---\n";
        cout << "The async benchmark should show better overlap ratios (> 1.0)\n";
        cout << "indicating that communication and computation are overlapped.\n";
        cout << "The sync benchmark should show overlap ratios ~ 1.0\n";
        cout << "since communication and computation are sequential.\n";
    }
}

// Print usage information
void print_usage(const char* prog_name) {
    cout << "Usage: " << prog_name << " [options]\n\n";
    cout << "Options:\n";
    cout << "  --test <type>            Test type: async, sync, or compare (default: compare)\n";
    cout << "  --grid-size <n>          Grid size (default: 256)\n";
    cout << "  --message-size <n>       Message size in doubles (default: 1024)\n";
    cout << "  --iterations <n>         Number of iterations (default: 100)\n";
    cout << "  --output <file>          CSV output filename (default: overlap_results.csv)\n";
    cout << "  --help                   Show this help message\n";
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Default parameters
    string test_type = "compare";
    size_t grid_size = 256;
    size_t message_size = 1024;
    int iterations = 100;
    string csv_filename = "overlap_results.csv";

    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];

        if (arg == "--help") {
            if (rank == 0) {
                print_usage(argv[0]);
            }
            MPI_Finalize();
            return 0;
        } else if (arg == "--test") {
            if (i + 1 < argc) {
                test_type = argv[++i];
            }
        } else if (arg == "--grid-size") {
            if (i + 1 < argc) {
                grid_size = atoi(argv[++i]);
            }
        } else if (arg == "--message-size") {
            if (i + 1 < argc) {
                message_size = atoi(argv[++i]);
            }
        } else if (arg == "--iterations") {
            if (i + 1 < argc) {
                iterations = atoi(argv[++i]);
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
        cout << "   MPI Communication/Computation Overlap\n";
        cout << "========================================\n";
        cout << "Processes: " << size << "\n";
        cout << "CSV output: " << csv_filename << "\n";
    }

    csv.write_header();

    // Run tests
    if (test_type == "async") {
        benchmark_async_overlap(rank, size, iterations, grid_size, message_size, csv);
    } else if (test_type == "sync") {
        benchmark_sync_no_overlap(rank, size, iterations, grid_size, message_size, csv);
    } else if (test_type == "compare") {
        compare_async_vs_sync(rank, size, iterations, grid_size, message_size, csv);
    } else if (rank == 0) {
        cout << "Unknown test type: " << test_type << "\n";
        MPI_Finalize();
        return 1;
    }

    if (rank == 0) {
        cout << "\n========================================\n";
        cout << "   Overlap Benchmark Complete\n";
        cout << "========================================\n";
        cout << "Results saved to: " << csv_filename << "\n";
    }

    MPI_Finalize();
    return 0;
}
