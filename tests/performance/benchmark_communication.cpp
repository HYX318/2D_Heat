/**
 * MPI Communication Performance Benchmark
 *
 * This file implements comprehensive benchmarks for MPI communication patterns:
 * - Ping-Pong (latency and bandwidth)
 * - Allreduce (collective reduction operations)
 * - Broadcast (tree-based broadcasting)
 * - Barrier (global synchronization)
 */

#include <mpi.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <algorithm>

using namespace std;

// Helper class for timing measurements
class Timer {
public:
    Timer() : start_time_(MPI_Wtime()) {}

    double elapsed() const {
        return MPI_Wtime() - start_time_;
    }

    void reset() {
        start_time_ = MPI_Wtime();
    }

private:
    double start_time_;
};

// Helper class for CSV output
class CSVWriter {
public:
    CSVWriter(int rank, const std::string& filename)
        : rank_(rank), filename_(filename) {}

    void write_header() {
        if (rank_ == 0) {
            ofstream file(filename_);
            if (file.is_open()) {
                file << "test_type,processes,message_size,time_us,bandwidth_mbps,operation\n";
                file.close();
            }
        }
    }

    void write(const std::string& test_type, int processes, size_t message_size,
               double time_us, double bandwidth_mbps, const std::string& operation = "") {
        if (rank_ == 0) {
            ofstream file(filename_, ios::app);
            if (file.is_open()) {
                file << test_type << "," << processes << "," << message_size << ","
                     << fixed << setprecision(4) << time_us << ","
                     << fixed << setprecision(4) << bandwidth_mbps << ","
                     << operation << "\n";
                file.close();
            }
        }
    }

private:
    int rank_;
    std::string filename_;
};

// Format bytes to human-readable string
string format_bytes(size_t bytes) {
    const char* units[] = {"B", "KB", "MB", "GB", "TB"};
    int unit_index = 0;
    double size = static_cast<double>(bytes);

    while (size >= 1024.0 && unit_index < 4) {
        size /= 1024.0;
        unit_index++;
    }

    ostringstream oss;
    oss << fixed << setprecision(1) << size << units[unit_index];
    return oss.str();
}

// Benchmark ping-pong communication using MPI_Send/MPI_Recv
void benchmark_pingpong_sendrecv(int rank, int size, int iterations,
                                   const vector<size_t>& message_sizes,
                                   CSVWriter& csv) {
    if (size < 2) {
        if (rank == 0) {
            cout << "\n=== Ping-Pong Benchmark (Send/Recv) ===\n";
            cout << "Error: Need at least 2 processes for ping-pong test\n";
        }
        return;
    }

    if (rank == 0) {
        cout << "\n=== Ping-Pong Benchmark (Send/Recv) ===\n";
        cout << "Iterations: " << iterations << "\n\n";
    }

    for (size_t msg_size : message_sizes) {
        // Allocate buffer
        vector<char> buffer(msg_size, static_cast<char>(rank));

        // Warm-up
        for (int i = 0; i < 10; i++) {
            if (rank == 0) {
                MPI_Send(buffer.data(), msg_size, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
                MPI_Recv(buffer.data(), msg_size, MPI_BYTE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else if (rank == 1) {
                MPI_Recv(buffer.data(), msg_size, MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(buffer.data(), msg_size, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        // Benchmark
        double total_time = 0.0;

        if (rank == 0 || rank == 1) {
            Timer timer;

            for (int i = 0; i < iterations; i++) {
                if (rank == 0) {
                    MPI_Send(buffer.data(), msg_size, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
                    MPI_Recv(buffer.data(), msg_size, MPI_BYTE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                } else {
                    MPI_Recv(buffer.data(), msg_size, MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(buffer.data(), msg_size, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
                }
            }

            total_time = timer.elapsed();
        }

        // Report results
        if (rank == 0) {
            double avg_time_us = (total_time / iterations) * 1e6;
            double time_per_message_us = avg_time_us / 2.0; // Round trip
            double bandwidth_mbps = (msg_size * 2.0) / (total_time / iterations) * 1e-6;

            cout << setw(10) << left << format_bytes(msg_size) << ": "
                 << setw(12) << right << fixed << setprecision(2) << time_per_message_us << " µs (latency)"
                 << setw(12) << right << fixed << setprecision(2) << bandwidth_mbps << " MB/s (bandwidth)\n";

            csv.write("pingpong_sendrecv", 2, msg_size, time_per_message_us, bandwidth_mbps, "round-trip");
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
}

// Benchmark ping-pong using MPI_Sendrecv
void benchmark_pingpong_sendrecv_combined(int rank, int size, int iterations,
                                            const vector<size_t>& message_sizes,
                                            CSVWriter& csv) {
    if (size < 2) {
        if (rank == 0) {
            cout << "\n=== Ping-Pong Benchmark (Sendrecv) ===\n";
            cout << "Error: Need at least 2 processes for ping-pong test\n";
        }
        return;
    }

    if (rank == 0) {
        cout << "\n=== Ping-Pong Benchmark (Sendrecv) ===\n";
        cout << "Iterations: " << iterations << "\n\n";
    }

    for (size_t msg_size : message_sizes) {
        vector<char> send_buffer(msg_size, static_cast<char>(rank));
        vector<char> recv_buffer(msg_size, 0);

        // Warm-up
        for (int i = 0; i < 10; i++) {
            int dest = (rank + 1) % 2;
            int source = (rank + 1) % 2;
            MPI_Sendrecv(send_buffer.data(), msg_size, MPI_BYTE, dest, 0,
                        recv_buffer.data(), msg_size, MPI_BYTE, source, 0,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        // Benchmark
        double total_time = 0.0;

        if (rank == 0 || rank == 1) {
            Timer timer;

            for (int i = 0; i < iterations; i++) {
                int dest = (rank + 1) % 2;
                int source = (rank + 1) % 2;
                MPI_Sendrecv(send_buffer.data(), msg_size, MPI_BYTE, dest, 0,
                            recv_buffer.data(), msg_size, MPI_BYTE, source, 0,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            total_time = timer.elapsed();
        }

        // Report results
        if (rank == 0) {
            double avg_time_us = (total_time / iterations) * 1e6;
            double bandwidth_mbps = (msg_size * 2.0) / (total_time / iterations) * 1e-6;

            cout << setw(10) << left << format_bytes(msg_size) << ": "
                 << setw(12) << right << fixed << setprecision(2) << avg_time_us << " µs"
                 << setw(12) << right << fixed << setprecision(2) << bandwidth_mbps << " MB/s\n";

            csv.write("pingpong_sendrecv_combined", 2, msg_size, avg_time_us, bandwidth_mbps, "sendrecv");
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
}

// Benchmark allreduce with different operations
void benchmark_allreduce(int rank, int size, int iterations,
                         const vector<size_t>& array_sizes,
                         CSVWriter& csv) {
    if (rank == 0) {
        cout << "\n=== Allreduce Benchmark ===\n";
        cout << "Iterations: " << iterations << "\n";
        cout << "Processes: " << size << "\n\n";
    }

    const char* op_names[] = {"SUM", "MAX", "MIN"};
    MPI_Op ops[] = {MPI_SUM, MPI_MAX, MPI_MIN};

    for (int op_idx = 0; op_idx < 3; op_idx++) {
        if (rank == 0) {
            cout << "\n--- Operation: " << op_names[op_idx] << " ---\n";
        }

        for (size_t array_size : array_sizes) {
            vector<double> send_buf(array_size, static_cast<double>(rank));
            vector<double> recv_buf(array_size, 0.0);

            // Warm-up
            for (int i = 0; i < 10; i++) {
                MPI_Allreduce(send_buf.data(), recv_buf.data(), array_size,
                             MPI_DOUBLE, ops[op_idx], MPI_COMM_WORLD);
            }

            MPI_Barrier(MPI_COMM_WORLD);

            // Benchmark
            Timer timer;

            for (int i = 0; i < iterations; i++) {
                MPI_Allreduce(send_buf.data(), recv_buf.data(), array_size,
                             MPI_DOUBLE, ops[op_idx], MPI_COMM_WORLD);
            }

            double total_time = timer.elapsed();
            double avg_time_us = (total_time / iterations) * 1e6;
            double bytes_per_iter = array_size * sizeof(double) * 2; // send + receive
            double bandwidth_mbps = bytes_per_iter / (total_time / iterations) * 1e-6;

            if (rank == 0) {
                cout << setw(10) << left << format_bytes(array_size * sizeof(double)) << ": "
                     << setw(12) << right << fixed << setprecision(2) << avg_time_us << " µs"
                     << setw(12) << right << fixed << setprecision(2) << bandwidth_mbps << " MB/s\n";

                csv.write("allreduce", size, array_size * sizeof(double), avg_time_us, bandwidth_mbps, op_names[op_idx]);
            }

            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
}

// Benchmark broadcast from root to all processes
void benchmark_broadcast(int rank, int size, int iterations,
                         const vector<size_t>& message_sizes,
                         CSVWriter& csv) {
    if (rank == 0) {
        cout << "\n=== Broadcast Benchmark ===\n";
        cout << "Iterations: " << iterations << "\n";
        cout << "Processes: " << size << "\n\n";
    }

    for (size_t msg_size : message_sizes) {
        vector<char> buffer(msg_size, static_cast<char>(rank));

        // Warm-up
        for (int i = 0; i < 10; i++) {
            MPI_Bcast(buffer.data(), msg_size, MPI_BYTE, 0, MPI_COMM_WORLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        // Benchmark
        Timer timer;

        for (int i = 0; i < iterations; i++) {
            MPI_Bcast(buffer.data(), msg_size, MPI_BYTE, 0, MPI_COMM_WORLD);
        }

        double total_time = timer.elapsed();
        double avg_time_us = (total_time / iterations) * 1e6;
        double bandwidth_mbps = msg_size / (total_time / iterations) * 1e-6;

        if (rank == 0) {
            cout << setw(10) << left << format_bytes(msg_size) << ": "
                 << setw(12) << right << fixed << setprecision(2) << avg_time_us << " µs"
                 << setw(12) << right << fixed << setprecision(2) << bandwidth_mbps << " MB/s\n";

            csv.write("broadcast", size, msg_size, avg_time_us, bandwidth_mbps, "from_root");
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
}

// Benchmark barrier (global synchronization)
void benchmark_barrier(int rank, int size, int iterations, CSVWriter& csv) {
    if (rank == 0) {
        cout << "\n=== Barrier Benchmark ===\n";
        cout << "Iterations: " << iterations << "\n";
        cout << "Processes: " << size << "\n\n";
    }

    // Warm-up
    for (int i = 0; i < 10; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Benchmark
    Timer timer;

    for (int i = 0; i < iterations; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
    }

    double total_time = timer.elapsed();
    double avg_time_us = (total_time / iterations) * 1e6;

    if (rank == 0) {
        cout << "Average barrier time: " << fixed << setprecision(2) << avg_time_us << " µs\n";
        cout << "Barrier frequency: " << fixed << setprecision(0) << (1.0 / (avg_time_us * 1e-6)) << " Hz\n";

        csv.write("barrier", size, 0, avg_time_us, 0.0, "synchronization");
    }
}

// Print usage information
void print_usage(const char* prog_name) {
    cout << "Usage: " << prog_name << " [options]\n\n";
    cout << "Options:\n";
    cout << "  --all                    Run all benchmarks\n";
    cout << "  --test <name>            Run specific test (pingpong, allreduce, broadcast, barrier)\n";
    cout << "  --iterations <n>         Number of iterations per test (default: 1000)\n";
    cout << "  --output <file>          CSV output filename (default: benchmark_results.csv)\n";
    cout << "  --help                   Show this help message\n";
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Default parameters
    int iterations = 1000;
    string csv_filename = "benchmark_results.csv";
    bool run_all = false;
    bool run_pingpong = false;
    bool run_allreduce = false;
    bool run_broadcast = false;
    bool run_barrier = false;

    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];

        if (arg == "--help") {
            if (rank == 0) {
                print_usage(argv[0]);
            }
            MPI_Finalize();
            return 0;
        } else if (arg == "--all") {
            run_all = true;
        } else if (arg == "--test") {
            if (i + 1 < argc) {
                string test_name = argv[++i];
                if (test_name == "pingpong") {
                    run_pingpong = true;
                } else if (test_name == "allreduce") {
                    run_allreduce = true;
                } else if (test_name == "broadcast") {
                    run_broadcast = true;
                } else if (test_name == "barrier") {
                    run_barrier = true;
                } else if (rank == 0) {
                    cout << "Unknown test: " << test_name << "\n";
                    MPI_Finalize();
                    return 1;
                }
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

    // If no specific tests specified, run all
    if (!run_all && !run_pingpong && !run_allreduce && !run_broadcast && !run_barrier) {
        run_all = true;
    }

    if (run_all) {
        run_pingpong = true;
        run_allreduce = true;
        run_broadcast = true;
        run_barrier = true;
    }

    // Test parameters
    vector<size_t> message_sizes = {
        1,                      // 1 byte
        1024,                   // 1 KB
        1024 * 1024,            // 1 MB
        10 * 1024 * 1024        // 10 MB
    };

    vector<size_t> array_sizes = {
        256,                    // 256 doubles (2 KB)
        1024,                   // 1024 doubles (8 KB)
        8192,                   // 8192 doubles (64 KB)
        65536,                  // 65536 doubles (512 KB)
        524288                  // 524288 doubles (4 MB)
    };

    CSVWriter csv(rank, csv_filename);

    if (rank == 0) {
        cout << "\n========================================\n";
        cout << "   MPI Communication Performance Benchmark\n";
        cout << "========================================\n";
        cout << "Number of processes: " << size << "\n";
        cout << "Iterations per test: " << iterations << "\n";
        cout << "CSV output: " << csv_filename << "\n";
    }

    csv.write_header();

    // Run benchmarks
    if (run_pingpong) {
        benchmark_pingpong_sendrecv(rank, size, iterations, message_sizes, csv);
        benchmark_pingpong_sendrecv_combined(rank, size, iterations, message_sizes, csv);
    }

    if (run_allreduce) {
        benchmark_allreduce(rank, size, iterations, array_sizes, csv);
    }

    if (run_broadcast) {
        benchmark_broadcast(rank, size, iterations, message_sizes, csv);
    }

    if (run_barrier) {
        benchmark_barrier(rank, size, iterations, csv);
    }

    if (rank == 0) {
        cout << "\n========================================\n";
        cout << "   Benchmark Complete\n";
        cout << "========================================\n";
        cout << "Results saved to: " << csv_filename << "\n";
    }

    MPI_Finalize();
    return 0;
}
