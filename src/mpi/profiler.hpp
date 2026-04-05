#ifndef PROFILER_HPP
#define PROFILER_HPP

#include <mpi.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <mutex>
#include <memory>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>

/**
 * @enum PhaseType
 * @brief Type of execution phase for timing
 */
enum class PhaseType {
    COMMUNICATION,  ///< Communication phase
    COMPUTATION,     ///< Computation phase
    IDLE           ///< Idle/waiting phase
};

/**
 * @struct CommunicationEvent
 * @brief Record of a single communication event
 */
struct CommunicationEvent {
    int source_or_dest;  ///< Source (for recv) or destination (for send)
    int bytes;           ///< Number of bytes transferred
    std::string tag;     ///< User-defined tag for the operation
    double start_time;   ///< Start time (MPI_Wtime)
    double end_time;     ///< End time (MPI_Wtime)
    bool is_send;        ///< true for send, false for receive

    CommunicationEvent(int src_dest, int b, const std::string& t,
                       double start, double end, bool send)
        : source_or_dest(src_dest), bytes(b), tag(t),
          start_time(start), end_time(end), is_send(send) {}
};

/**
 * @struct TimingRegion
 * @brief Timing information for a named region
 */
struct TimingRegion {
    double total_time;        ///< Total time spent in this region
    double start_time;        ///< Start time of current region (if active)
    size_t call_count;        ///< Number of times region was entered
    bool is_active;           ///< Whether region is currently active

    TimingRegion() : total_time(0.0), start_time(0.0),
                     call_count(0), is_active(false) {}
};

/**
 * @struct CommunicationStats
 * @brief Statistics for a specific communication type
 */
struct CommunicationStats {
    size_t send_count;        ///< Number of sends
    size_t recv_count;        ///< Number of receives
    size_t total_bytes_sent;  ///< Total bytes sent
    size_t total_bytes_recv;  ///< Total bytes received
    double total_time;        ///< Total time spent in this type

    CommunicationStats() : send_count(0), recv_count(0),
                          total_bytes_sent(0), total_bytes_recv(0),
                          total_time(0.0) {}
};

/**
 * @class Profiler
 * @brief MPI performance profiling and analysis tool
 *
 * This class provides comprehensive profiling for MPI applications, including:
 * - Hierarchical timing (total, communication, computation, idle)
 * - Communication tracking (message counts, data volumes, operation types)
 * - Performance metrics (communication/computation ratio, load balance)
 * - Detailed reporting (console, CSV, JSON)
 *
 * Features:
 * - Zero-overhead when disabled (compile-time option)
 * - Thread-safe statistics
 * - MPI Profiling Interface (PMPI) support
 * - Integration with MPI_Wtime() for accurate timing
 *
 * Usage example:
 * @code
 * Profiler profiler(true);  // Enable PMPI profiling
 * profiler.start_computation("compute_step");
 * // ... do computation ...
 * profiler.end_computation("compute_step");
 * profiler.start_communication("exchange_ghost");
 * MPI_Isend(...);
 * MPI_Irecv(...);
 * MPI_Waitall(...);
 * profiler.end_communication("exchange_ghost");
 * profiler.print_report();
 * @endcode
 */
class Profiler {
public:
    /**
     * @brief Construct Profiler
     * @param enable_mpi_pmpi Enable MPI Profiling Interface (P MPI)
     * @param mpi_rank Rank of this process (for reporting)
     * @param mpi_size Total number of processes
     *
     * If enable_mpi_pmpi is true, attempts to enable MPI profiling.
     */
    explicit Profiler(bool enable_mpi_pmpi = true,
                      int mpi_rank = 0,
                      int mpi_size = 1);

    /**
     * @brief Destructor - automatically write reports if configured
     */
    ~Profiler();

    // Disable copying
    Profiler(const Profiler&) = delete;
    Profiler& operator=(const Profiler&) = delete;

    // Allow moving
    Profiler(Profiler&&) noexcept = default;
    Profiler& operator=(Profiler&&) noexcept = default;

    /**
     * @brief Start timing a communication region
     * @param name Name/identifier for the region
     *
     * Starts timing a named communication phase. The region must be ended
     * with a matching call to end_communication().
     */
    void start_communication(const std::string& name);

    /**
     * @brief End timing a communication region
     * @param name Name/identifier for the region (must match start)
     *
     * Ends timing a named communication phase and updates statistics.
     */
    void end_communication(const std::string& name);

    /**
     * *brief Start timing a computation region
     * @param name Name/identifier for the region
     *
     * Starts timing a named computation phase. The region must be ended
     * with a matching call to end_computation().
     */
    void start_computation(const std::string& name);

    /**
     * @brief End timing a computation region
     * @param name Name/identifier for the region (must match start)
     *
     * Ends timing a named computation phase and updates statistics.
     */
    void end_computation(const std::string& name);

    /**
     * @brief Record a send operation
     * @param dest Destination rank
     * @param bytes Number of bytes sent
     * @param tag User-defined tag for the operation
     *
     * Records statistics about a send operation. Call this after each MPI send.
     */
    void record_send(int dest, int bytes, const std::string& tag);

    /**
     * @brief Record a receive operation
     * @param src Source rank
     * @param bytes Number of bytes received
     * @param tag User-defined tag for the operation
     *
     * Records statistics about a receive operation. Call this after each MPI receive.
     */
    void record_recv(int src, int bytes, const std::string& tag);

    /**
     * @brief Print a human-readable performance report to stdout
     *
     * Prints detailed timing and communication statistics.
     */
    void print_report() const;

    /**
     * @brief Write performance report to a file
     * @param filename Output filename
     *
     * Writes a human-readable report to the specified file.
     */
    void write_report(const std::string& filename) const;

    /**
     * @brief Write CSV format report
     * @param filename Output filename
     *
     * Writes timing and communication data in CSV format for analysis.
     */
    void write_csv(const std::string& filename) const;

    /**
     * @brief Write JSON format report
     * @param filename Output filename
     *
     * Writes timing and communication data in JSON format.
     */
    void write_json(const std::string& filename) const;

    // Accessors for statistics

    /**
     * @brief Get total elapsed time since profiler creation
     * @return Total time in seconds
     */
    double total_time() const { return MPI_Wtime() - start_time_; }

    /**
     * @brief Get total communication time
     * @return Communication time in seconds
     */
    double communication_time() const { return total_comm_time_; }

    /**
     * @brief Get total computation time
     * @return Computation time in seconds
     */
    double computation_time() const { return total_comp_time_; }

    /**
     * @brief Get total idle time
     * @return Idle time in seconds
     */
    double idle_time() const { return total_idle_time_; }

    /**
     * @brief Get total number of messages (sends + receives)
     * @return Total message count
     */
    size_t total_messages() const {
        return global_stats_.send_count + global_stats_.recv_count;
    }

    /**
     * @brief Get total bytes transferred (sent + received)
     * @return Total bytes
     */
    size_t total_bytes() const {
        return global_stats_.total_bytes_sent + global_stats_.total_bytes_recv;
    }

    /**
     * @brief Get communication/computation ratio
     * @return Ratio (0 if no computation)
     */
    double comm_comp_ratio() const {
        return total_comp_time_ > 0 ? total_comm_time_ / total_comp_time_ : 0.0;
    }

    /**
     * @brief Get percentage of time spent in communication
     * @return Percentage (0-100)
     */
    double comm_percentage() const {
        double total = total_time();
        return total > 0 ? (total_comm_time_ / total) * 100.0 : 0.0;
    }

    /**
     * @brief Get percentage of time spent in computation
     * @return Percentage (0-100)
     */
    double comp_percentage() const {
        double total = total_time();
        return total > 0 ? (total_comp_time_ / total) * 100.0 : 0.0;
    }

    /**
     * @brief Get timing statistics for a specific region
     * @param name Region name
     * @return Const pointer to TimingRegion, or nullptr if not found
     */
    const TimingRegion* get_region_stats(const std::string& name) const;

    /**
     * @brief Get communication statistics for a specific tag
     * @param tag Communication operation tag
     * @return Const pointer to CommunicationStats, or nullptr if not found
     */
    const CommunicationStats* get_comm_stats(const std::string& tag) const;

    /**
     * @brief Get all recorded communication events
     * @return Const reference to event list
     */
    const std::vector<CommunicationEvent>& events() const { return events_; }

    /**
     * @brief Enable or disable profiling
     * @param enabled true to enable, false to disable
     *
     * When disabled, all profiling operations become no-ops (zero overhead).
     */
    void set_enabled(bool enabled) { enabled_ = enabled; }

    /**
     * @brief Check if profiling is enabled
     * @return true if enabled, false otherwise
     */
    bool is_enabled() const { return enabled_; }

    /**
     * @brief Reset all statistics
     *
     * Clears all timing regions, communication statistics, and events.
     */
    void reset();

    /**
     * @brief Synchronize timing across all processes
     *
     * Computes load balance statistics based on timing from all processes.
     * Requires MPI to be initialized.
     */
    void synchronize_and_analyze();

    /**
     * @brief Get load balance efficiency
     * @return Efficiency (0-1), or -1 if not computed
     */
    double load_balance_efficiency() const { return load_balance_efficiency_; }

    /**
     * @brief Get parallel efficiency
     * @param serial_time Time for serial execution
     * @param num_procs Number of processors
     * @return Efficiency (0-1), or -1 if not computed
     */
    double parallel_efficiency(double serial_time, int num_procs) const;

private:
    bool enabled_;                    ///< Whether profiling is enabled
    bool mpi_pmpi_enabled_;           ///< Whether PMPI profiling is enabled
    int mpi_rank_;                    ///< MPI rank of this process
    int mpi_size_;                    ///< Total number of MPI processes
    double start_time_;               ///< Profiler start time

    // Timing statistics
    double total_comm_time_;          ///< Total communication time
    double total_comp_time_;          ///< Total computation time
    double total_idle_time_;          ///< Total idle time
    double current_phase_start_;      ///< Start time of current phase
    PhaseType current_phase_;         ///< Current phase type

    // Region timing
    std::unordered_map<std::string, TimingRegion> timing_regions_;

    // Communication statistics
    CommunicationStats global_stats_;      ///< Global communication stats
    std::unordered_map<std::string, CommunicationStats> comm_stats_;  ///< Per-tag stats

    // Event tracking
    std::vector<CommunicationEvent> events;  ///< List of all communication events
    std::vector<CommunicationEvent> events_; ///< Alias for compatibility

    // Load balance
    double load_balance_efficiency_;         ///< Load balance efficiency (0-1)

    // Thread safety (optional, for MPI_THREAD_MULTIPLE)
    mutable std::mutex mutex_;               ///< Mutex for thread-safe updates

    /**
     * @brief Lock mutex for thread-safe operations
     * @return Lock guard
     */
    std::lock_guard<std::mutex> lock() const { return std::lock_guard<std::mutex>(mutex_); }

    /**
     * @brief Get current time using MPI_Wtime
     * @return Time in seconds
     */
    double get_time() const { return MPI_Wtime(); }

    /**
     * @brief Format time with appropriate units
     * @param seconds Time in seconds
     * @return Formatted string
     */
    static std::string format_time(double seconds);

    /**
     * @brief Format bytes with appropriate units
     * @param bytes Number of bytes
     * @return Formatted string
     */
    static std::string format_bytes(size_t bytes);

    /**
     * @brief Write header for CSV file
     * @param file Output file stream
     */
    void write_csv_header(std::ofstream& file) const;

    /**
     * @brief Write timing data to CSV file
     * @param file Output file stream
     */
    void write_csv_timing(std::ofstream& file) const;

    /**
     * @brief Write communication data to CSV file
     * @param file Output file stream
     */
    void write_csv_communication(std::ofstream& file) const;
};

#endif // PROFILER_HPP
