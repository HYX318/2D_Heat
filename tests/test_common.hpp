/**
 * @file test_common.hpp
 * @brief Common test utilities and macros for Google Test framework
 */

#ifndef TEST_COMMON_HPP
#define TEST_COMMON_HPP

#include <gtest/gtest.h>
#include <mpi.h>
#include <memory>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <functional>

/**
 * @namespace test_utils
 * @brief Utility functions and classes for testing
 */
namespace test_utils {

/**
 * @brief Check if MPI is initialized
 * @return true if MPI is initialized, false otherwise
 */
inline bool is_mpi_initialized() {
    int initialized = 0;
    MPI_Initialized(&initialized);
    return initialized != 0;
}

/**
 * @brief Check if MPI is finalized
 * @return true if MPI is finalized, false otherwise
 */
inline bool is_mpi_finalized() {
    int finalized = 0;
    MPI_Finalized(&finalized);
    return finalized != 0;
}

/**
 * @brief Initialize MPI if not already initialized
 * @param argc Argument count from main()
 * @param argv Argument vector from main()
 */
inline void ensure_mpi_initialized(int& argc, char** argv) {
    if (!is_mpi_initialized()) {
        MPI_Init(&argc, &argv);
    }
}

/**
 * @brief Get current MPI rank
 * @return Current process rank, or -1 if MPI not initialized
 */
inline int get_mpi_rank() {
    if (!is_mpi_initialized()) {
        return -1;
    }
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

/**
 * @brief Get total number of MPI processes
 * @return Number of processes, or -1 if MPI not initialized
 */
inline int get_mpi_size() {
    if (!is_mpi_initialized()) {
        return -1;
    }
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

/**
 * @brief Check if current process is root (rank 0)
 * @return true if rank == 0, false otherwise
 */
inline bool is_root() {
    return get_mpi_rank() == 0;
}

/**
 * @brief MPI-safe print that only outputs from root process
 * @param message Message to print
 */
inline void mpi_print(const std::string& message) {
    if (is_root()) {
        std::cout << message << std::endl;
    }
}

/**
 * @brief Helper class to measure test execution time
 */
class Timer {
public:
    Timer() : start_(std::chrono::high_resolution_clock::now()) {}

    /**
     * @brief Get elapsed time in milliseconds
     * @return Elapsed time in milliseconds
     */
    double elapsed_ms() const {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::milli>(end - start_).count();
    }

    /**
     * @brief Get elapsed time in seconds
     * @return Elapsed time in seconds
     */
    double elapsed_s() const {
        return elapsed_ms() / 1000.0;
    }

    /**
     * @brief Reset the timer
     */
    void reset() {
        start_ = std::chrono::high_resolution_clock::now();
    }

private:
    std::chrono::high_resolution_clock::time_point start_;
};

/**
 * @brief RAII helper for MPI initialization/finalization in tests
 */
class MPIGuard {
public:
    MPIGuard(int& argc, char** argv, bool finalize_on_destruct = false)
        : should_finalize_(finalize_on_destruct) {
        was_initialized_ = is_mpi_initialized();
        if (!was_initialized_) {
            MPI_Init(&argc, &argv);
        }
    }

    ~MPIGuard() {
        if (should_finalize_ && !was_initialized_ && !is_mpi_finalized()) {
            MPI_Finalize();
        }
    }

    // Prevent copying
    MPIGuard(const MPIGuard&) = delete;
    MPIGuard& operator=(const MPIGuard&) = delete;

private:
    bool was_initialized_;
    bool should_finalize_;
};

/**
 * @brief Skip test if MPI is not initialized
 */
inline void skip_if_no_mpi() {
    if (!is_mpi_initialized()) {
        GTEST_SKIP() << "MPI not initialized. Run with mpirun/mpiexec.";
    }
}

/**
 * @brief Skip test if MPI is not initialized or process count doesn't match
 * @param required_procs Required number of processes
 */
inline void skip_if_mpi_procs_not(int required_procs) {
    skip_if_no_mpi();

    int size = get_mpi_size();
    if (size != required_procs) {
        std::stringstream ss;
        ss << "This test requires exactly " << required_procs
           << " MPI processes, but running with " << size << " processes.";
        GTEST_SKIP() << ss.str();
    }
}

/**
 * @brief Skip test if MPI is not initialized or minimum process count not met
 * @param min_procs Minimum number of processes required
 */
inline void skip_if_mpi_procs_less_than(int min_procs) {
    skip_if_no_mpi();

    int size = get_mpi_size();
    if (size < min_procs) {
        std::stringstream ss;
        ss << "This test requires at least " << min_procs
           << " MPI processes, but running with " << size << " processes.";
        GTEST_SKIP() << ss.str();
    }
}

/**
 * @brief Helper to execute a function only on root process
 * @param func Function to execute on root
 */
inline void execute_on_root(std::function<void()> func) {
    if (is_root()) {
        func();
    }
}

/**
 * @brief Synchronize all processes with a barrier
 */
inline void mpi_barrier() {
    if (is_mpi_initialized() && !is_mpi_finalized()) {
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

/**
 * @brief Format double with specified precision
 * @param value Value to format
 * @param precision Number of decimal places
 * @return Formatted string
 */
inline std::string format_double(double value, int precision = 6) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(precision) << value;
    return ss.str();
}

/**
 * @brief Compare two doubles with relative tolerance
 * @param a First value
 * @param b Second value
 * @param rel_tol Relative tolerance
 * @return true if values are approximately equal
 */
inline bool approx_equal(double a, double b, double rel_tol = 1e-6) {
    return std::abs(a - b) <= rel_tol * std::max(std::abs(a), std::abs(b));
}

} // namespace test_utils

#define skip_if_no_mpi() \
    do { \
        if (!test_utils::is_mpi_initialized()) { \
            GTEST_SKIP() << "MPI not initialized. Run with mpirun/mpiexec."; \
        } \
    } while(0)

#define skip_if_mpi_procs_not(required_procs) \
    do { \
        skip_if_no_mpi(); \
        int mpi_size = test_utils::get_mpi_size(); \
        if (mpi_size != (required_procs)) { \
            std::stringstream ss; \
            ss << "This test requires exactly " << (required_procs) \
               << " MPI processes, but running with " << mpi_size << " processes."; \
            GTEST_SKIP() << ss.str(); \
        } \
    } while(0)

#define skip_if_mpi_procs_less_than(min_procs) \
    do { \
        skip_if_no_mpi(); \
        int mpi_size = test_utils::get_mpi_size(); \
        if (mpi_size < (min_procs)) { \
            std::stringstream ss; \
            ss << "This test requires at least " << (min_procs) \
               << " MPI processes, but running with " << mpi_size << " processes."; \
            GTEST_SKIP() << ss.str(); \
        } \
    } while(0)

/**
 * @brief Custom assertion macro for MPI initialization check
 */
#define ASSERT_MPI_INITIALIZED() \
    do { \
        ASSERT_TRUE(test_utils::is_mpi_initialized()) \
            << "MPI must be initialized for this test"; \
    } while(0)

/**
 *   @brief Custom assertion macro for MPI process count
 */
#define ASSERT_MPI_PROCS(required) \
    do { \
        ASSERT_MPI_INITIALIZED(); \
        int size = test_utils::get_mpi_size(); \
        ASSERT_EQ(size, required) \
            << "Test requires " << required << " processes, got " << size; \
    } while(0)

/**
 * @brief Custom assertion macro for checking if process is root
 */
#define ASSERT_IS_ROOT() \
    do { \
        ASSERT_MPI_INITIALIZED(); \
        ASSERT_EQ(test_utils::get_mpi_rank(), 0) \
            << "This assertion must be run on root process"; \
    } while(0)

/**
 * @brief Macro to execute code only on root process
 */
#define ON_ROOT(code) \
    do { \
        if (test_utils::is_root()) { \
            code; \
        } \
    } while(0)

/**
 * @brief Macro to measure and print test execution time
 */
#define MEASURE_TIME() \
    test_utils::Timer _test_timer

/**
 * @brief Macro to print elapsed time
 */
#define PRINT_ELAPSED_TIME() \
    std::cout << "Test execution time: " \
              << test_utils::format_double(_test_timer.elapsed_ms(), 3) \
              << " ms" << std::endl

#endif // TEST_COMMON_HPP
