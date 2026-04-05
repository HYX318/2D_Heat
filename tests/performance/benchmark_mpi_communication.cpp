#include <mpi.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <cstring>

/**
 * @file benchmark_mpi_communication.cpp
 * @brief Comprehensive MPI communication performance benchmarks
 *
 * This program runs various MPI communication benchmarks to measure:
 * - Ghost cell exchange performance
 * - Reduction operation performance
 * - Broadcast performance
 * - All-to-all communication performance
 *
 * Results are written in multiple formats:
 * - CSV files for plotting and analysis
 * - JSON files for programmatic processing
 * - Console summary for quick review
 */

struct BenchmarkResult {
    std::string name;
    int num_procs;
    int message_size;
    int iterations;
    double total_time;
    double avg_time;
    double min_time;
    double max_time;
    double throughput_mb_s;
    double latency_ns;

    BenchmarkResult() : num_procs(0), message_size(0), iterations(0),
                       total_time(0.0), avg_time(0.0), min_time(0.0),
                       max_time(0.0), throughput_mb_s(0.0), latency_ns(0.0) {}
};

class Stopwatch {
public:
    void start() {
        start_time_ = MPI_Wtime();
    }

    void stop() {
        elapsed_ = MPI_Wtime() - start_time_;
    }

    double elapsed() const {
        return elapsed_;
    }

    void reset() {
        elapsed_ = 0.0;
    }

private:
    double start_time_;
    double elapsed_;
};

/**
 * @brief Benchmark blocking send/recv for ghost cell exchange
 */
BenchmarkResult benchmark_ghost_cell_blocking(int rank, int size,
                                               int message_size, int iterations) {
    BenchmarkResult result;
    result.name = "ghost_cell_blocking";
    result.num_procs = size;
    result.message_size = message_size;
    result.iterations = iterations;

    std::vector<double> data(message_size);
    std::fill(data.begin(), data.end(), rank + 1.0);

    std::vector<double> recv_buffer(message_size);
    int partner = (rank + 1) % size;

    // Warmup
    if (rank < size / 2) {
        MPI_Send(data.data(), message_size, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD);
    } else {
        MPI_Recv(recv_buffer.data(), message_size, MPI_DOUBLE,
                 partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    std::vector<double> timings;

    for (int iter = 0; iter < iterations; ++iter) {
        MPI_Barrier(MPI_COMM_WORLD);

        Stopwatch sw;
        sw.start();

        if (rank < size / 2) {
            MPI_Send(data.data(), message_size, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD);
        } else {
            MPI_Recv(recv_buffer.data(), message_size, MPI_DOUBLE,
                     partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        sw.stop();

        timings.push_back(sw.elapsed());
    }

    result.total_time = std::accumulate(timings.begin(), timings.end(), 0.0);
    result.avg_time = result.total_time / iterations;
    result.min_time = *std::min_element(timings.begin(), timings.end());
    result.max_time = *std::max_element(timings.begin(), timings.end());

    // Compute latency (ping-pong style)
    result.latency_ns = (result.avg_time * 1e9) / 2.0;

    // Compute throughput (for large messages)
    double bytes_per_sec = (message_size * sizeof(double)) / result.avg_time;
    result.throughput_mb_s = bytes_per_sec / (1024.0 * 1024.0);

    return result;
}

/**
 * @brief Benchmark non-blocking send/recv (asynchronous)
 */
BenchmarkResult benchmark_ghost_cell_async(int rank, int size,
                                             int message_size, int iterations) {
    BenchmarkResult result;
    result.name = "ghost_cell_async";
    result.num_procs = size;
    result.message_size = message_size;
    result.iterations = iterations;

    std::vector<double> data(message_size);
    std::fill(data.begin(), data.end(), rank + 1.0);

    std::vector<double> recv_buffer(message_size);
    int partner = (rank + 1) % size;

    std::vector<double> timings;

    for (int iter = 0; iter < iterations; ++iter) {
        MPI_Barrier(MPI_COMM_WORLD);

        Stopwatch sw;
        sw.start();

        MPI_Request request;
        MPI_Status status;

        if (rank < size / 2) {
            MPI_Isend(data.data(), message_size, MPI_DOUBLE, partner, 0,
                      MPI_COMM_WORLD, &request);
        } else {
            MPI_Irecv(recv_buffer.data(), message_size, MPI_DOUBLE,
                      partner, 0, MPI_COMM_WORLD, &request);
        }

        MPI_Wait(&request, &status);

        sw.stop();

        timings.push_back(sw.elapsed());
    }

    result.total_time = std::accumulate(timings.begin(), timings.end(), 0.0);
    result.avg_time = result.total_time / iterations;
    result.min_time = *std::min_element(timings.begin(), timings.end());
    result.max_time = *std::max_element(timings.begin(), timings.end());

    result.latency_ns = (result.avg_time * 1e9) / 2.0;
    double bytes_per_sec = (message_size * sizeof(double)) / result.avg_time;
    result.throughput_mb_s = bytes_per_sec / (1024.0 * 1024.0);

    return result;
}

/**
 * @brief Benchmark MPI_Sendrecv for simultaneous exchange
 */
BenchmarkResult benchmark_sendrecv(int rank, int size,
                                     int message_size, int iterations) {
    BenchmarkResult result;
    result.name = "sendrecv";
    result.num_procs = size;
    result.message_size = message_size;
    result.iterations = iterations;

    std::vector<double> send_data(message_size);
    std::fill(send_data.begin(), send_data.end(), rank + 1.0);

    std::vector<double> recv_data(message_size);
    int dest = (rank + 1) % size;
    int source = (rank + size - 1) % size;

    std::vector<double> timings;

    for (int iter = 0; iter < iterations; ++iter) {
        MPI_Barrier(MPI_COMM_WORLD);

        Stopwatch sw;
        sw.start();

        MPI_Sendrecv(send_data.data(), message_size, MPI_DOUBLE, dest, 0,
                     recv_data.data(), message_size, MPI_DOUBLE, source, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        sw.stop();

        timings.push_back(sw.elapsed());
    }

    result.total_time = std::accumulate(timings.begin(), timings.end(), 0.0);
    result.avg_time = result.total_time / iterations;
    result.min_time = *std::min_element(timings.begin(), timings.end());
    result.max_time = *std::max_element(timings.begin(), timings.end());

    result.latency_ns = result.avg_time * 1e9;
    double bytes_per_sec = (2 * message_size * sizeof(double)) / result.avg_time;
    result.throughput_mb_s = bytes_per_sec / (1024.0 * 1024.0);

    return result;
}

/**
 * @brief Benchmark MPI_Reduce
 */
BenchmarkResult benchmark_reduce(int rank, int size,
                                  int message_size, int iterations) {
    BenchmarkResult result;
    result.name = "reduce_sum";
    result.num_procs = size;
    result.message_size = message_size;
    result.iterations = iterations;

    std::vector<double> data(message_size);
    std::fill(data.begin(), data.end(), rank + 1.0);

    std::vector<double> result_data(message_size);

    std::vector<double> timings;

    for (int iter = 0; iter < iterations; ++iter) {
        MPI_Barrier(MPI_COMM_WORLD);

        Stopwatch sw;
        sw.start();

        MPI_Reduce(data.data(), result_data.data(), message_size,
                   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        sw.stop();

        timings.push_back(sw.elapsed());
    }

    result.total_time = std::accumulate(timings.begin(), timings.end(), 0.0);
    result.avg_time = result.total_time / iterations;
    result.min_time = *std::min_element(timings.begin(), timings.end());
    result.max_time = *std::max_element(timings.begin(), timings.end());

    // Reduce is a collective operation - latency doesn't apply directly
    result.latency_ns = result.avg_time * 1e9;

    // Throughput for reduction (each process contributes data)
    double total_bytes = message_size * sizeof(double) * size;
    double bytes_per_sec = total_bytes / result.avg_time;
    result.throughput_mb_s = bytes_per_sec / (1024.0 * 1024.0);

    return result;
}

/**
 * @brief Benchmark MPI_Allreduce
 */
BenchmarkResult benchmark_allreduce(int rank, int size,
                                     int message_size, int iterations) {
    BenchmarkResult result;
    result.name = "allreduce_sum";
    result.num_procs = size;
    result.message_size = message_size;
    result.iterations = iterations;

    std::vector<double> data(message_size);
    std::fill(data.begin(), data.end(), rank + 1.0);

    std::vector<double> result_data(message_size);

    std::vector<double> timings;

    for (int iter = 0; iter < iterations; ++iter) {
        MPI_Barrier(MPI_COMM_WORLD);

        Stopwatch sw;
        sw.start();

        MPI_Allreduce(data.data(), result_data.data(), message_size,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        sw.stop();

        timings.push_back(sw.elapsed());
    }

    result.total_time = std::accumulate(timings.begin(), timings.end(), 0.0);
    result.avg_time = result.total_time / iterations;
    result.min_time = *std::min_element(timings.begin(), timings.end());
    result.max_time = *std::max_element(timings.begin(), timings.end());

    result.latency_ns = result.avg_time * 1e9;
    double total_bytes = message_size * sizeof(double) * size;
    double bytes_per_sec = total_bytes / result.avg_time;
    result.throughput_mb_s = bytes_per_sec / (1024.0 * 1024.0);

    return result;
}

/**
 * @brief Benchmark MPI_Bcast
 */
BenchmarkResult benchmark_bcast(int rank, int size,
                                 int message_size, int iterations) {
    BenchmarkResult result;
    result.name = "broadcast";
    result.num_procs = size;
    result.message_size = message_size;
    result.iterations = iterations;

    std::vector<double> data(message_size);

    if (rank == 0) {
        std::fill(data.begin(), data.end(), 1.0);
    }

    std::vector<double> timings;

    for (int iter = 0; iter < iterations; ++iter) {
        MPI_Barrier(MPI_COMM_WORLD);

        Stopwatch sw;
        sw.start();

        MPI_Bcast(data.data(), message_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        sw.stop();

        timings.push_back(sw.elapsed());
    }

    result.total_time = std::accumulate(timings.begin(), timings.end(), 0.0);
    result.avg_time = result.total_time / iterations;
    result.min_time = *std::min_element(timings.begin(), timings.end());
    result.max_time = *std::max_element(timings.begin(), timings.end());

    // For broadcast, latency is approximately measured
    result.latency_ns = result.avg_time * 1e9;

    // Throughput: data size / time
    double bytes_per_sec = (message_size * sizeof(double)) / result.avg_time;
    result.throughput_mb_s = bytes_per_sec / (1024.0 * 1024.0);

    return result;
}

/**
 * @brief Benchmark MPI_Alltoall
 */
BenchmarkResult benchmark_alltoall(int rank, int size,
                                    int message_size, int iterations) {
    BenchmarkResult result;
    result.name = "alltoall";
    result.num_procs = size;
    result.message_size = message_size;
    result.iterations = iterations;

    int send_count = message_size;
    int recv_count = message_size;

    std::vector<double> send_data(send_count * size);
    std::vector<double> recv_data(recv_count * size);

    // Initialize send data
    for (int i = 0; i < size; ++i) {
        std::fill(send_data.begin() + i * send_count,
                   send_data.begin() + (i + 1) * send_count,
                   rank + 1.0);
    }

    std::vector<double> timings;

    for (int iter = 0; iter < iterations; ++iter) {
        MPI_Barrier(MPI_COMM_WORLD);

        Stopwatch sw;
        sw.start();

        MPI_Alltoall(send_data.data(), send_count, MPI_DOUBLE,
                     recv_data.data(), recv_count, MPI_DOUBLE,
                     MPI_COMM_WORLD);

        sw.stop();

        timings.push_back(sw.elapsed());
    }

    result.total_time = std::accumulate(timings.begin(), timings.end(), 0.0);
    result.avg_time = result.total_time / iterations;
    result.min_time = *std::min_element(timings.begin(), timings.end());
    result.max_time = *std::max_element(timings.begin(), timings.end());

    result.latency_ns = result.avg_time * 1e9;

    // Total data transferred per process
    double total_bytes = send_count * recv_count * size * sizeof(double);
    double bytes_per_sec = total_bytes / result.avg_time;
    result.throughput_mb_s = bytes_per_sec / (1024.0 * 1024.0);

    return result;
}

/**
 * @brief Print benchmark result to console
 */
void print_benchmark_result(const BenchmarkResult& result, int rank) {
    if (rank == 0) {
        std::cout << std::left << std::setw(20) << result.name
                  << std::setw(10) << result.message_size
                  << std::setw(10) << result.iterations
                  << std::setw(12) << std::fixed << std::setprecision(6) << result.avg_time
                  << std::setw(12) << std::setprecision(6) << result.min_time
                  << std::setw(12) << std::setprecision(6) << result.max_time
                  << std::setw(15) << std::setprecision(2) << result.throughput_mb_s
                  << std::setw(12) << std::setprecision(1) << result.latency_ns
                  << "\n";
    }
}

/**
 * @brief Write benchmark result to CSV file
 */
void write_csv_header(std::ofstream& file) {
    file << "benchmark_name,num_processes,message_size,iterations,"
         << "avg_time_s,min_time_s,max_time_s,throughput_mb_s,latency_ns\n";
}

void write_csv_result(std::ofstream& file, const BenchmarkResult& result) {
    file << result.name << ","
         << result.num_procs << ","
         << result.message_size << ","
         << result.iterations << ","
         << result.avg_time << ","
         << result.min_time << ","
         << result.max_time << ","
         << result.throughput_mb_s << ","
         << result.latency_ns << "\n";
}

/**
 * @brief Write benchmark result to JSON file
 */
void write_json_header(std::ofstream& file, int rank, int size) {
    file << "{\n";
    file << "  \"processes\": " << size << ",\n";
    file << "  \"results\": [\n";
}

void write_json_result(std::ofstream& file, const BenchmarkResult& result,
                       bool is_last) {
    file << "    {\n";
    file << "      \"name\": \"" << result.name << "\",\n";
    file << "      \"message_size\": " << result.message_size << ",\n";
    file << "      \"iterations\": " << result.iterations << ",\n";
    file << "      \"avg_time_s\": " << result.avg_time << ",\n";
    file << "      \"min_time_s\": " << result.min_time << ",\n";
    file << "      \"max_time_s\": " << result.max_time << ",\n";
    file << "      \"throughput_mb_s\": " << result.throughput_mb_s << ",\n";
    file << "      \"latency_ns\": " << result.latency_ns << "\n";
    file << "    }" << (is_last ? "\n" : ",\n");
}

void write_json_footer(std::ofstream& file) {
    file << "  ]\n";
    file << "}\n";
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Parse command line arguments
    int warmup_iters = 10;
    int benchmark_iters = 100;
    std::vector<int> message_sizes = {1, 10, 100, 1000, 10000};
    std::string output_prefix = "benchmark_results";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--iters" && i + 1 < argc) {
            benchmark_iters = std::stoi(argv[++i]);
        } else if (arg == "--sizes" && i + 1 < argc) {
            message_sizes.clear();
            std::istringstream iss(argv[++i]);
            std::string size_str;
            while (std::getline(iss, size_str, ',')) {
                message_sizes.push_back(std::stoi(size_str));
            }
        } else if (arg == "--output" && i + 1 < argc) {
            output_prefix = argv[++i];
        }
    }

    if (rank == 0) {
        std::cout << "\n===========================================\n";
        std::cout << "MPI Communication Benchmarks\n";
        std::cout << "===========================================\n";
        std::cout << "Processes: " << size << "\n";
        std::cout << "Iterations: " << benchmark_iters << "\n";
        std::cout << "Message sizes: ";
        for (size_t i = 0; i < message_sizes.size(); ++i) {
            std::cout << message_sizes[i];
            if (i < message_sizes.size() - 1) std::cout << ", ";
        }
        std::cout << "\n";
        std::cout << "===========================================\n\n";
    }

    std::vector<BenchmarkResult> all_results;

    // Print header
    if (rank == 0) {
        std::cout << std::left << std::setw(20) << "Benchmark"
                  << std::setw(10) << "Msg Size"
                  << std::setw(10) << "Iters"
                  << std::setw(12) << "Avg (s)"
                  << std::setw(12) << "Min (s)"
                  << std::setw(12) << "Max (s)"
                  << std::setw(15) << "Throughput"
                  << std::setw(12) << "Latency"
                  << "\n";
        std::cout << std::string(100, '-') << "\n";
    }

    // Run benchmarks for each message size
    for (int msg_size : message_sizes) {
        // Ghost cell exchange benchmarks
        if (size >= 2) {
            all_results.push_back(
                benchmark_ghost_cell_blocking(rank, size, msg_size, benchmark_iters));
            print_benchmark_result(all_results.back(), rank);

            all_results.push_back(
                benchmark_ghost_cell_async(rank, size, msg_size, benchmark_iters));
            print_benchmark_result(all_results.back(), rank);

            all_results.push_back(
                benchmark_sendrecv(rank, size, msg_size, benchmark_iters));
            print_benchmark_result(all_results.back(), rank);
        }

        // Reduction benchmarks
        all_results.push_back(
            benchmark_reduce(rank, size, msg_size, benchmark_iters));
        print_benchmark_result(all_results.back(), rank);

        all_results.push_back(
            benchmark_allreduce(rank, size, msg_size, benchmark_iters));
        print_benchmark_result(all_results.back(), rank);

        // Broadcast benchmark
        all_results.push_back(
            benchmark_bcast(rank, size, msg_size, benchmark_iters));
        print_benchmark_result(all_results.back(), rank);

        // All-to-all benchmark
        all_results.push_back(
            benchmark_alltoall(rank, size, msg_size, benchmark_iters));
        print_benchmark_result(all_results.back(), rank);

        if (rank == 0) {
            std::cout << "\n";
        }
    }

    // Write CSV results
    if (rank == 0) {
        std::string csv_filename = output_prefix + ".csv";
        std::ofstream csv_file(csv_filename);
        if (csv_file) {
            write_csv_header(csv_file);
            for (const auto& result : all_results) {
                write_csv_result(csv_file, result);
            }
            csv_file.close();
            std::cout << "CSV results written to " << csv_filename << "\n";
        }

        // Write JSON results
        std::string json_filename = output_prefix + ".json";
        std::ofstream json_file(json_filename);
        if (json_file) {
            write_json_header(json_file, rank, size);
            for (size_t i = 0; i < all_results.size(); ++i) {
                write_json_result(json_file, all_results[i],
                                  i == all_results.size() - 1);
            }
            write_json_footer(json_file);
            json_file.close();
            std::cout << "JSON results written to " << json_filename << "\n";
        }
    }

    if (rank == 0) {
        std::cout << "\n===========================================\n";
        std::cout << "Benchmarking Complete\n";
        std::cout << "===========================================\n\n";
    }

    MPI_Finalize();
    return 0;
}
