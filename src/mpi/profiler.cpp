#include "profiler.hpp"
#include <iostream>
#include <cmath>
#include <iomanip>

Profiler::Profiler(bool enable_mpi_pmpi, int mpi_rank, int mpi_size)
    : enabled_(true)
    , mpi_pmpi_enabled_(false)
    , mpi_rank_(mpi_rank)
    , mpi_size_(mpi_size)
    , start_time_(MPI_Wtime())
    , total_comm_time_(0.0)
    , total_comp_time_(0.0)
    , total_idle_time_(0.0)
    , current_phase_start_(0.0)
    , current_phase_(PhaseType::IDLE)
    , load_balance_efficiency_(-1.0)
{
    if (enable_mpi_pmpi) {
        // Attempt to enable MPI profiling interface
        // Note: MPI_Pcontrol might not be available in all MPI implementations
        // This is a best-effort attempt
        int flag = 0;
        MPI_Initialized(&flag);
        if (flag) {
            // Try to enable profiling (implementation-specific)
            // MPI_Pcontrol(1);  // Uncomment if your MPI supports this
            mpi_pmpi_enabled_ = true;
        }
    }
}

Profiler::~Profiler() {
    // Automatically write report if needed
    // User can call write_report() explicitly with their preferred filename
}

void Profiler::start_communication(const std::string& name) {
    if (!enabled_) return;

    auto lk = lock();
    double now = get_time();

    // End current phase if active
    if (current_phase_ != PhaseType::IDLE) {
        double elapsed = now - current_phase_start_;
        if (current_phase_ == PhaseType::COMPUTATION) {
            total_comp_time_ += elapsed;
        } else if (current_phase_ == PhaseType::COMMUNICATION) {
            total_comm_time_ += elapsed;
        }
    }

    // Start new communication phase
    current_phase_ = PhaseType::COMMUNICATION;
    current_phase_start_ = now;

    // Start timing region
    auto& region = timing_regions_[name];
    if (!region.is_active) {
        region.start_time = now;
        region.is_active = true;
    }
}

void Profiler::end_communication(const std::string& name) {
    if (!enabled_) return;

    auto lk = lock();
    double now = get_time();

    // End current communication phase
    if (current_phase_ == PhaseType::COMMUNICATION) {
        double elapsed = now - current_phase_start_;
        total_comm_time_ += elapsed;
        current_phase_ = PhaseType::IDLE;
    }

    // End timing region
    auto it = timing_regions_.find(name);
    if (it != timing_regions_.end() && it->second.is_active) {
        double elapsed = now - it->second.start_time;
        it->second.total_time += elapsed;
        it->second.call_count++;
        it->second.is_active = false;
    }
}

void Profiler::start_computation(const std::string& name) {
    if (!enabled_) return;

    auto lk = lock();
    double now = get_time();

    // End current phase if active
    if (current_phase_ != PhaseType::IDLE) {
        double elapsed = now - current_phase_start_;
        if (current_phase_ == PhaseType::COMPUTATION) {
            total_comp_time_ += elapsed;
        } else if (current_phase_ == PhaseType::COMMUNICATION) {
            total_comm_time_ += elapsed;
        }
    }

    // Start new computation phase
    current_phase_ = PhaseType::COMPUTATION;
    current_phase_start_ = now;

    // Start timing region
    auto& region = timing_regions_[name];
    if (!region.is_active) {
        region.start_time = now;
        region.is_active = true;
    }
}

void Profiler::end_computation(const std::string& name) {
    if (!enabled_) return;

    auto lk = lock();
    double now = get_time();

    // End current computation phase
    if (current_phase_ == PhaseType::COMPUTATION) {
        double elapsed = now - current_phase_start_;
        total_comp_time_ += elapsed;
        current_phase_ = PhaseType::IDLE;
    }

    // End timing region
    auto it = timing_regions_.find(name);
    if (it != timing_regions_.end() && it->second.is_active) {
        double elapsed = now - it->second.start_time;
        it->second.total_time += elapsed;
        it->second.call_count++;
        it->second.is_active = false;
    }
}

void Profiler::record_send(int dest, int bytes, const std::string& tag) {
    if (!enabled_) return;

    auto lk = lock();
    double now = get_time();

    // Update global stats
    global_stats_.send_count++;
    global_stats_.total_bytes_sent += bytes;

    // Update per-tag stats
    auto& stats = comm_stats_[tag];
    stats.send_count++;
    stats.total_bytes_sent += bytes;

    // Record event
    CommunicationEvent event(dest, bytes, tag, now, now, true);
    events_.push_back(event);
}

void Profiler::record_recv(int src, int bytes, const std::string& tag) {
    if (!enabled_) return;

    auto lk = lock();
    double now = get_time();

    // Update global stats
    global_stats_.recv_count++;
    global_stats_.total_bytes_recv += bytes;

    // Update per-tag stats
    auto& stats = comm_stats_[tag];
    stats.recv_count++;
    stats.total_bytes_recv += bytes;

    // Record event
    CommunicationEvent event(src, bytes, tag, now, now, false);
    events_.push_back(event);
}

void Profiler::print_report() const {
    if (!enabled_) {
        std::cout << "Profiling is disabled.\n";
        return;
    }

    std::cout << "\n===========================================\n";
    std::cout << "MPI Performance Profiler Report\n";
    std::cout << "===========================================\n";
    std::cout << "Process: " << mpi_rank_ << " / " << mpi_size_ << "\n";
    std::cout << "-------------------------------------------\n\n";

    // Overall timing
    std::cout << "Overall Timing:\n";
    std::cout << "  Total Time:       " << format_time(total_time()) << "\n";
    std::cout << "  Communication:    " << format_time(total_comm_time_)
              << " (" << std::fixed << std::setprecision(2) << comm_percentage() << "%)\n";
    std::cout << "  Computation:      " << format_time(total_comp_time_)
              << " (" << comp_percentage() << "%)\n";
    std::cout << "  Comm/Comp Ratio:  " << std::fixed << std::setprecision(3)
              << comm_comp_ratio() << "\n";

    if (load_balance_efficiency_ >= 0.0) {
        std::cout << "  Load Balance:     " << std::fixed << std::setprecision(1)
                  << (load_balance_efficiency_ * 100.0) << "%\n";
    }
    std::cout << "\n";

    // Communication statistics
    std::cout << "Communication Statistics:\n";
    std::cout << "  Total Messages:   " << total_messages() << "\n";
    std::cout << "  Sends:           " << global_stats_.send_count << "\n";
    std::cout << "  Receives:        " << global_stats_.recv_count << "\n";
    std::cout << "  Total Bytes:     " << format_bytes(total_bytes()) << "\n";
    std::cout << "  Bytes Sent:      " << format_bytes(global_stats_.total_bytes_sent) << "\n";
    std::cout << "  Bytes Received:  " << format_bytes(global_stats_.total_bytes_recv) << "\n";

    if (global_stats_.send_count > 0) {
        double avg_send = static_cast<double>(global_stats_.total_bytes_sent) / global_stats_.send_count;
        std::cout << "  Avg Send Size:   " << format_bytes(static_cast<size_t>(avg_send)) << "\n";
    }
    if (global_stats_.recv_count > 0) {
        double avg_recv = static_cast<double>(global_stats_.total_bytes_recv) / global_stats_.recv_count;
        std::cout << "  Avg Recv Size:   " << format_bytes(static_cast<size_t>(avg_recv)) << "\n";
    }
    std::cout << "\n";

    // Timing regions
    if (!timing_regions_.empty()) {
        std::cout << "Timing Regions:\n";
        std::cout << std::left << std::setw(30) << "  Region"
                  << std::right << std::setw(12) << "Calls"
                  << std::setw(15) << "Total Time"
                  << std::setw(15) << "Avg/Call"
                  << "\n";
        std::cout << std::string(72, '-') << "\n";

        for (const auto& pair : timing_regions_) {
            const std::string& name = pair.first;
            const TimingRegion& region = pair.second;
            double avg = region.call_count > 0 ? region.total_time / region.call_count : 0.0;

            std::cout << "  " << std::left << std::setw(28) << name
                      << std::right << std::setw(12) << region.call_count
                      << std::setw(15) << format_time(region.total_time)
                      << std::setw(15) << format_time(avg)
                      << "\n";
        }
        std::cout << "\n";
    }

    // Per-tag communication statistics
    if (!comm_stats_.empty()) {
        std::cout << "Communication by Tag:\n";
        std::cout << std::left << std::setw(30) << "  Tag"
                  << std::right << std::setw(10) << "Sends"
                  << std::setw(10) << "Recvs"
                  << std::setw(15) << "Bytes Sent"
                  << std::setw(15) << "Bytes Recv"
                  << "\n";
        std::cout << std::string(80, '-') << "\n";

        for (const auto& pair : comm_stats_) {
            const std::string& tag = pair.first;
            const CommunicationStats& stats = pair.second;

            std::cout << "  " << std::left << std::setw(28) << tag
                      << std::right << std::setw(10) << stats.send_count
                      << std::setw(10) << stats.recv_count
                      << std::setw(15) << format_bytes(stats.total_bytes_sent)
                      << std::setw(15) << format_bytes(stats.total_bytes_recv)
                      << "\n";
        }
        std::cout << "\n";
    }

    std::cout << "===========================================\n\n";
}

void Profiler::write_report() const {
    // Use default filename if called without argument
    std::string filename = "profiler_report_rank_" + std::to_string(mpi_rank_) + ".txt";
    write_report(filename);
}

void Profiler::write_report(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    // Redirect cout to file for print_report
    std::streambuf* old_cout = std::cout.rdbuf();
    std::cout.rdbuf(file.rdbuf());

    print_report();

    std::cout.rdbuf(old_cout);
    file.close();

    std::cout << "Report written to " << filename << "\n";
}

void Profiler::write_csv(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    file << "# MPI Performance Profiler CSV Report\n";
    file << "# Process: " << mpi_rank_ << " / " << mpi_size_ << "\n";
    file << "# Total Time: " << total_time() << "\n";
    file << "\n";

    // Write timing regions
    write_csv_header(file);
    write_csv_timing(file);
    file << "\n";

    // Write communication data
    file << "# Communication by Tag\n";
    file << "Tag,Sends,Receives,Bytes_Sent,Bytes_Received\n";
    for (const auto& pair : comm_stats_) {
        const std::string& tag = pair.first;
        const CommunicationStats& stats = pair.second;
        file << "\"" << tag << "\","
             << stats.send_count << ","
             << stats.recv_count << ","
             << stats.total_bytes_sent << ","
             << stats.total_bytes_recv << "\n";
    }

    file.close();
    std::cout << "CSV report written to " << filename << "\n";
}

void Profiler::write_csv_header(std::ofstream& file) const {
    file << "# Timing Regions\n";
    file << "Region,Calls,Total_Time,Avg_Time_Per_Call\n";
}

void Profiler::write_csv_timing(std::ofstream& file) const {
    for (const auto& pair : timing_regions_) {
        const std::string& name = pair.first;
        const TimingRegion& region = pair.second;
        double avg = region.call_count > 0 ? region.total_time / region.call_count : 0.0;

        file << "\"" << name << "\","
             << region.call_count << ","
             << region.total_time << ","
             << avg << "\n";
    }
}

void Profiler::write_csv_communication(std::ofstream& file) const {
    file << "# Communication Events\n";
    file << "Timestamp,Source_Dest,Bytes,Tag,Is_Send\n";
    for (const auto& event : events_) {
        file << event.start_time << ","
             << event.source_or_dest << ","
             << event.bytes << ","
             << "\"" << event.tag << "\","
             << (event.is_send ? "true" : "false") << "\n";
    }
}

void Profiler::write_json(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    file << "{\n";
    file << "  \"process_rank\": " << mpi_rank_ << ",\n";
    file << "  \"num_processes\": " << mpi_size_ << ",\n";
    file << "  \"total_time\": " << total_time() << ",\n";

    // Timing
    file << "  \"timing\": {\n";
    file << "    \"communication_time\": " << total_comm_time_ << ",\n";
    file << "    \"computation_time\": " << total_comp_time_ << ",\n";
    file << "    \"idle_time\": " << total_idle_time_ << ",\n";
    file << "    \"comm_percentage\": " << comm_percentage() << ",\n";
    file << "    \"comp_percentage\": " << comp_percentage() << ",\n";
    file << "    \"comm_comp_ratio\": " << comm_comp_ratio() << "\n";
    file << "  },\n";

    // Communication stats
    file << "  \"communication\": {\n";
    file << "    \"total_messages\": " << total_messages() << ",\n";
    file << "    \"send_count\": " << global_stats_.send_count << ",\n";
    file << "    \"recv_count\": " << global_stats_.recv_count << ",\n";
    file << "    \"total_bytes\": " << total_bytes() << ",\n";
    file << "    \"bytes_sent\": " << global_stats_.total_bytes_sent << ",\n";
    file << "    \"bytes_received\": " << global_stats_.total_bytes_recv << "\n";
    file << "  },\n";

    // Timing regions
    file << "  \"timing_regions\": {\n";
    bool first = true;
    for (const auto& pair : timing_regions_) {
        const std::string& name = pair.first;
        const TimingRegion& region = pair.second;
        double avg = region.call_count > 0 ? region.total_time / region.call_count : 0.0;

        if (!first) file << ",\n";
        first = false;

        file << "    \"" << name << "\": {\n";
        file << "      \"total_time\": " << region.total_time << ",\n";
        file << "      \"call_count\": " << region.call_count << ",\n";
        file << "      \"avg_time\": " << avg << "\n";
        file << "    }";
    }
    file << "\n  },\n";

    // Per-tag stats
    file << "  \"communication_by_tag\": {\n";
    first = true;
    for (const auto& pair : comm_stats_) {
        const std::string& tag = pair.first;
        const CommunicationStats& stats = pair.second;

        if (!first) file << ",\n";
        first = false;

        file << "    \"" << tag << "\": {\n";
        file << "      \"send_count\": " << stats.send_count << ",\n";
        file << "      "recv_count\": " << stats.recv_count << ",\n";
        file << "      \"bytes_sent\": " << stats.total_bytes_sent << ",\n";
        file << "      \"bytes_received\": " << stats.total_bytes_recv << "\n";
        file << "    }";
    }
    file << "\n  }\n";

    file << "}\n";
    file.close();

    std::cout << "JSON report written to " << filename << "\n";
}

const TimingRegion* Profiler::get_region_stats(const std::string& name) const {
    auto it = timing_regions_.find(name);
    return (it != timing_regions_.end()) ? &it->second : nullptr;
}

const CommunicationStats* Profiler::get_comm_stats(const std::string& tag) const {
    auto it = comm_stats_.find(tag);
    return (it != comm_stats_.end()) ? &it->second : nullptr;
}

void Profiler::reset() {
    auto lk = lock();
    start_time_ = get_time();
    total_comm_time_ = 0.0;
    total_comp_time_ = 0.0;
    total_idle_time_ = 0.0;
    current_phase_start_ = 0.0;
    current_phase_ = PhaseType::IDLE;
    timing_regions_.clear();
    global_stats_ = CommunicationStats();
    comm_stats_.clear();
    events_.clear();
    load_balance_efficiency_ = -1.0;
}

void Profiler::synchronize_and_analyze() {
    if (mpi_size_ <= 1) {
        load_balance_efficiency_ = 1.0;
        return;
    }

    int flag = 0;
    MPI_Initialized(&flag);
    if (!flag) return;

    // Compute local total time
    double local_total = total_comm_time();

    // Gather all times
    std::vector<double> all_times(mpi_size_);
    MPI_Gather(&local_total, 1, MPI_DOUBLE,
               all_times.data(), 1, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    // Root process computes load balance efficiency
    if (mpi_rank_ == 0) {
        double max_time = *std::max_element(all_times.begin(), all_times.end());
        double avg_time = std::accumulate(all_times.begin(), all_times.end(), 0.0) / mpi_size_;
        load_balance_efficiency_ = avg_time / max_time;
    }

    // Broadcast efficiency to all processes
    MPI_Bcast(&load_balance_efficiency_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

double Profiler::parallel_efficiency(double serial_time, int num_procs) const {
    if (num_procs <= 0 || serial_time <= 0) return -1.0;

    double parallel_time = total_time();
    double speedup = serial_time / parallel_time;
    double efficiency = speedup / num_procs;

    return std::min(efficiency, 1.0);
}

std::string Profiler::format_time(double seconds) {
    std::ostringstream oss;

    if (seconds < 1e-6) {
        oss << std::scientific << std::setprecision(2) << seconds * 1e9 << " ns";
    } else if (seconds < 1e-3) {
        oss << std::fixed << std::setprecision(2) << seconds * 1e6 << " us";
    } else if (seconds < 1.0) {
        oss << std::fixed << std::setprecision(2) << seconds * 1e3 << " ms";
    } else if (seconds < 60.0) {
        oss << std::fixed << std::setprecision(3) << seconds << " s";
    } else {
        int minutes = static_cast<int>(seconds / 60.0);
        double secs = seconds - minutes * 60.0;
        oss << std::fixed << std::setprecision(1) << minutes << "m " << secs << "s";
    }

    return oss.str();
}

std::string Profiler::format_bytes(size_t bytes) {
    std::ostringstream oss;

    if (bytes < 1024) {
        oss << bytes << " B";
    } else if (bytes < 1024 * 1024) {
        oss << std::fixed << std::setprecision(2) << (bytes / 1024.0) << " KB";
    } else if (bytes < 1024ull * 1024 * 1024) {
        oss << std::fixed << std::setprecision(2) << (bytes / (1024.0 * 1024.0)) << " MB";
    } else {
        oss << std::fixed << std::setprecision(2) << (bytes / (1024.0 * 1024.0 * 1024.0)) << " GB";
    }

    return oss.str();
}
