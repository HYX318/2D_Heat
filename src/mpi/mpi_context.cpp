#include "mpi_context.hpp"
#include <iostream>
#include <sstream>

MPIContext::MPIContext(int& argc, char** argv)
    : rank_(0), size_(1), initialized_by_us_(false), finalized_(false) {
    // Check if MPI is already initialized
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);

    if (!mpi_initialized) {
        // Request thread support - MPI_THREAD_FUNNELED is sufficient for most cases
        int provided_thread_support = 0;
        int mpi_init_result = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED,
                                               &provided_thread_support);

        if (mpi_init_result != MPI_SUCCESS) {
            std::ostringstream oss;
            oss << "MPI_Init_thread failed with error code: " << mpi_init_result;
            throw std::runtime_error(oss.str());
        }

        initialized_by_us_ = true;
    }

    // Get process rank and size
    int mpi_comm_rank_result = MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    if (mpi_comm_rank_result != MPI_SUCCESS) {
        std::ostringstream oss;
        oss << "MPI_Comm_rank failed with error code: " << mpi_comm_rank_result;
        throw std::runtime_error(oss.str());
    }

    int mpi_comm_size_result = MPI_Comm_size(MPI_COMM_WORLD, &size_);
    if (mpi_comm_size_result != MPI_SUCCESS) {
        std::ostringstream oss;
        oss << "MPI_Comm_size failed with error code: " << mpi_comm_size_result;
        throw std::runtime_error(oss.str());
    }
}

MPIContext::~MPIContext() {
    try {
        // Check if MPI is still initialized
        int mpi_initialized = 0;
        MPI_Initialized(&mpi_initialized);

        // Check if MPI is already finalized
        int mpi_finalized = 0;
        MPI_Finalized(&mpi_finalized);

        // Only finalize if we initialized it and it's not already finalized
        if (mpi_initialized && !mpi_finalized && initialized_by_us_ && !finalized_) {
            int mpi_finalize_result = MPI_Finalize();
            if (mpi_finalize_result != MPI_SUCCESS) {
                std::cerr << "Warning: MPI_Finalize failed with error code: "
                          << mpi_finalize_result << std::endl;
            }
            finalized_ = true;
        }
    } catch (...) {
        // Destructor should never throw - catch any exceptions
        std::cerr << "Warning: Exception caught in MPIContext destructor" << std::endl;
    }
}

MPIContext::MPIContext(MPIContext&& other) noexcept
    : rank_(other.rank_),
      size_(other.size_),
      initialized_by_us_(other.initialized_by_us_),
      finalized_(other.finalized_) {
    // Mark the moved-from object as invalid
    other.initialized_by_us_ = false;
    other.finalized_ = true;
    other.rank_ = -1;
    other.size_ = -1;
}

MPIContext& MPIContext::operator=(MPIContext&& other) noexcept {
    if (this != &other) {
        // Clean up existing resources if we own them
        if (initialized_by_us_ && !finalized_) {
            try {
                int mpi_initialized = 0;
                MPI_Initialized(&mpi_initialized);

                int mpi_finalized = 0;
                MPI_Finalized(&mpi_finalized);

                if (mpi_initialized && !mpi_finalized) {
                    MPI_Finalize();
                }
            } catch (...) {
                // Assignment operator should never throw
                std::cerr << "Warning: Exception caught in MPIContext move assignment"
                          << std::endl;
            }
        }

        // Transfer ownership
        rank_ = other.rank_;
        size_ = other.size_;
        initialized_by_us_ = other.initialized_by_us_;
        finalized_ = other.finalized_;

        // Mark the moved as invalid
        other.initialized_by_us_ = false;
        other.finalized_ = true;
        other.rank_ = -1;
        other.size_ = -1;
    }
    return *this;
}

void MPIContext::barrier() const {
    if (finalized_) {
        throw std::runtime_error("Cannot call barrier on finalized MPI context");
    }

    int mpi_barrier_result = MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_barrier_result != MPI_SUCCESS) {
        std::ostringstream oss;
        oss << "MPI_Barrier failed with error code: " << mpi_barrier_result;
        throw std::runtime_error(oss.str());
    }
}

void MPIContext::abort(int error_code, const std::string& message) const {
    if (is_root()) {
        std::cerr << "MPI Abort: " << message << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, error_code);
}
