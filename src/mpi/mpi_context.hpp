#ifndef MPI_CONTEXT_HPP
#define MPI_CONTEXT_HPP

#include <mpi.h>
#include <string>
#include <stdexcept>

/**
 * @class MPIContext
 * @brief RAII wrapper for MPI initialization and finalization
 *
 * This class provides a safe, exception-safe wrapper around MPI initialization
 * and finalization. It automatically handles MPI_Init and MPI_Finalize calls,
 * preventing common errors like forgetting to finalize MPI or multiple initialization.
 *
 * Features:
 * - RAII pattern: automatic cleanup on destruction
 * - Thread-safe: supports MPI thread initialization
 * - Move-only: prevents accidental copies, allows moves
 * - Exception-safe: handles errors gracefully
 */
class MPIContext {
public:
    /**
     * @brief Construct and initialize MPI
     * @param argc Argument count from main()
     * @param argv Argument vector from main()
     * @throws std::runtime_error if MPI initialization fails
     *
     * The constructor checks if MPI is already initialized. If not, it
     * initializes MPI with MPI_THREAD_FUNNELED support, retrieves the
     * process rank and size, and sets up the context.
     */
    MPIContext(int& argc, char** argv);

    /**
     * @brief Destructor - automatically finalize MPI
     *
     * Checks if MPI is still initialized and if this context is responsible
     * for initialization. If so, calls MPI_Finalize() to clean up.
     * Exception-safe: catches and logs any errors during finalization.
     */
    ~MPIContext();

    // Delete copy constructor and copy assignment
    MPIContext(const MPIContext&) = delete;
    MPIContext& operator=(const MPIContext&) = delete;

    /**
     * @brief Move constructor
     * @param other The MPIContext to move from
     *
     * Transfers ownership of MPI initialization from another context.
     * The moved-from context becomes invalid and should not be used.
     */
    MPIContext(MPIContext&& other) noexcept;

    /**
     * @brief Move assignment operator
     * @param other The MPIContext to move from
     * @return Reference to this object
     *
     * Transfers ownership of MPI initialization from another context.
     * If this context currently owns MPI, it is finalized first.
     * The moved-from context becomes invalid and should not be used.
     */
    MPIContext& operator=(MPIContext&& other) noexcept;

    /**
     * @brief Get the rank of the current process
     * @return Process rank (0 to size-1)
     */
    int rank() const { return rank_; }

    /**
     * @brief Get the total number of processes
     * @return Total number of processes in the MPI communicator
     */
    int size() const { return size_; }

    /**
     * @brief Check if this process is the root process
     * @return true if rank == 0, false otherwise
     */
    bool is_root() const { return rank_ == 0; }

    /**
     * @brief Block until all processes reach this point
     *
     * Calls MPI_Barrier to synchronize all processes in the communicator.
     * This is useful for coordinating timing or ensuring all processes
     * have completed a phase of computation.
     *
     * @throws std::runtime_error if the barrier call fails
     */
    void barrier() const;

    /**
     * @brief Abort the MPI execution
     * @param error_code Error code to return to the environment
     * @param message Error message to display
     *
     * Calls MPI_Abort to terminate all MPI processes.
     * This is typically used when a fatal error occurs that cannot be recovered from.
     */
    void abort(int error_code, const std::string& message) const;

private:
    int rank_;               ///< Process rank in MPI_COMM_WORLD
    int size_;               ///< Total number of processes in MPI_COMM_WORLD
    bool initialized_by_us_; ///< True if this context initialized MPI
    bool finalized_;         ///< True if MPI has been finalized
};

#endif // MPI_CONTEXT_HPP
