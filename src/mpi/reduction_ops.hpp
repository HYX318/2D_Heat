#ifndef REDUCTION_OPS_HPP
#define REDUCTION_OPS_HPP

#include <mpi.h>
#include <vector>
#include <type_traits>
#include <stdexcept>
#include <cstdint>
#include <limits>

/**
 * @file reduction_ops.hpp
 * @brief Type-safe MPI reduction operations
 *
 * This header defines the ReductionOps class which provides a type-safe
 * interface to MPI reduction operations, wrapping MPI_Reduce, MPI_Allreduce,
 * MPI_Scan, MPI_Exscan, and related operations.
 */

// Forward declaration to avoid circular dependency
class MPIContext;

/**
 * @enum ReductionOp
 * @brief Supported reduction operation types
 *
 * These map to standard MPI reduction operations, providing a type-safe
 * enumeration for reduction operations.
 */
enum class ReductionOp {
    SUM,      ///< Summation (MPI_SUM)
    MAX,      ///< Maximum (MPI_MAX)
    MIN,      ///< Minimum (MPI_MIN)
    PROD,     ///< Product (MPI_PROD)
    LAND,     ///< Logical AND (MPI_LAND)
    LOR,      ///< Logical OR (MPI_LOR)
    LXOR,     ///< Logical XOR (MPI_LXOR)
    BAND,     ///< Bitwise AND (MPI_BAND)
    BOR,      ///< Bitwise OR (MPI_BOR)
    BXOR      ///< Bitwise XOR (MPI_BXOR)
};

/**
 * @struct MaxLocResult
 * @brief Result of maxloc reduction operation
 *
 * Contains both the maximum value and the rank of the process
 * where it was found.
 */
template<typename T>
struct MaxLocResult {
    T value;      ///< The maximum value found
    int rank;     ///< Rank of the process containing the maximum
};

/**
 * @struct MinLocResult
 * @brief Result of minloc reduction operation
 *
 * Contains both the minimum value and the rank of the process
 * where it was found.
 */
template<typename T>
struct MinLocResult {
    T value;      ///< The minimum value found
    int rank;     ///< Rank of the process containing the minimum
};

/**
 * @class ReductionOps
 * @brief Type-safe wrapper for MPI reduction operations
 *
 * This class provides a convenient, type-safe interface to MPI reduction
 * operations. It supports both scalar and array reductions, as well as
 * special operations like scan, exscan, maxloc, minloc, and broadcast.
 *
 * Example usage:
 * @code
 * ReductionOps ops;  // Use MPI_COMM_WORLD
 *
 * // Simple reduction
 * int result = ops.reduce(my_value, ReductionOp::SUM, 0);
 *
 * // Allreduce (broadcast to all)
 * double sum = ops.allreduce(my_double, ReductionOp::SUM);
 *
 * // Array reduction
 * std::vector<double> result_array = ops.allreduce_array(values, ReductionOp::MAX);
 *
 * // Find max and its location
 * auto max_result = ops.maxloc(my_value);
 * @endcode
 *
 * Thread safety: This class is thread-safe as long as underlying MPI
 * implementation supports MPI_THREAD_MULTIPLE or MPI_THREAD_SERIALIZED.
 */
class ReductionOps {
public:
    /**
     * @brief Construct with default communicator (MPI_COMM_WORLD)
     *
     * Creates a ReductionOps object that operates on MPI_COMM_WORLD.
     */
    ReductionOps();

    /**
     * @brief Construct with custom communicator
     * @param comm MPI communicator to use for operations
     *
     * Creates a ReductionOps object that operates on the specified communicator.
     */
    explicit ReductionOps(MPI_Comm comm);

    /**
     * @brief Construct from MPIContext
     * @param ctx MPIContext providing the communicator
     *
     * Creates a ReductionOps object using the communicator from MPIContext.
     */
    explicit ReductionOps(const MPIContext& ctx);

    /**
     * @brief Destructor
     */
    ~ReductionOps() = default;

    // Delete copy operations (communicator cannot be copied)
    ReductionOps(const ReductionOps&) = delete;
    ReductionOps& operator=(const ReductionOps&) = delete;

    // Allow move operations
    ReductionOps(ReductionOps&&) noexcept = default;
    ReductionOps& operator=(ReductionOps&&) noexcept = default;

    /**
     * @brief Reduce scalar value to root process
     * @tparam T Type of value (must be: int, long, long long, unsigned, unsigned long, unsigned long long, float, double, long double)
     * @param value Value to reduce from this process
     * @param op Reduction operation to apply
     * @param root Rank of root process where result is stored
     * @return Reduced value on root process, undefined on other processes
     * @throws std::runtime_error if MPI operation fails
     * @throws std::invalid_argument if type T is not supported
     *
     * Performs a reduction operation across all processes in the communicator,
     * with the result available only on the root process.
     */
    template<typename T>
    T reduce(T value, ReductionOp op, int root) const;

    /**
     * @brief Allreduce scalar value (broadcast to all)
     * @tparam T Type of value (must be arithmetic type)
     * @param value Value to reduce from this process
     * @param op Reduction operation to apply
     * @return Reduced value (same on all processes)
     * @throws std::runtime_error if MPI operation fails
     * @throws std::invalid_argument if type T is not supported
     *
     * Performs a reduction operation across all processes and broadcasts
     * the result to all processes.
     */
    template<typename T>
    T allreduce(T value, ReductionOp op) const;

    /**
     * @brief Scan operation (inclusive prefix reduction)
     * @tparam T Type of value (must be arithmetic type)
     * @param value Value from this process
     * @param op Reduction operation to apply
     * @return Result of reduction on processes 0 through current rank
     * @throws std::runtime_error if MPI operation fails
     * @throws std::invalid_argument if type T is not supported
     *
     * Performs an inclusive scan, where each process receives the reduction
     * of values from processes 0 through its own rank.
     */
    template<typename T>
    T scan(T value, ReductionOp op) const;

    /**
     * @brief Exclusive scan operation
     * @tparam T Type of value (must be arithmetic type)
     * @param value Value from this process
     * @param op Reduction operation to apply
     * @return Result of reduction on processes 0 through (current rank - 1)
     * @throws std::runtime_error if MPI operation fails
     * @throws std::invalid_argument if type T is not supported
     *
     * Performs an exclusive scan, where each process receives the reduction
     * of values from processes 0 through (rank - 1). Process 0 receives an
     * identity value.
     */
    template<typename T>
    T exscan(T value, ReductionOp op) const;

    /**
     * @brief Reduce array to root process
     * @tparam T Type of array elements (must be arithmetic type)
     * @param values Array to reduce from this process
     * @param op Reduction operation to apply
     * @param root Rank of root process where result is stored
     * @return Reduced array on root process, empty on other processes
     * @throws std::runtime_error if MPI operation fails
     * @throws std::invalid_argument if type T is not supported or array sizes don't match
     *
     * Performs element-wise reduction across all processes. All processes
     * must provide arrays of the same size.
     */
    template<typename T>
    std::vector<T> reduce_array(const std::vector<T>& values,
                                ReductionOp op, int root) const;

    /**
     * @brief Allreduce array (broadcast to all)
     * @tparam T Type of array elements (must be arithmetic type)
     * @param values Array to reduce from this process
     * @param op Reduction operation to apply
     * @return Reduced array (same on all processes)
     * @throws std::runtime_error if MPI operation fails
     * @throws std::invalid_argument if type T is not supported or array sizes don't match
     *
     * Performs element-wise reduction across all processes and broadcasts
     * the result to all processes.
     */
    template<typename T>
    std::vector<T> allreduce_array(const std::vector<T>& values,
                                   ReductionOp op) const;

    /**
     * @brief Find maximum value and its location
     * @tparam T Type of value (must be arithmetic type)
     * @param value Value from this process
     * @return MaxLocResult containing max value and rank where it was found
     * @throws std::runtime_error if MPI operation fails
     * @throws std::invalid_argument if type T is not supported
     *
     * Finds the maximum value across all processes and returns both the
     * value and the rank of the process where it was found.
     */
    template<typename T>
    MaxLocResult<T> maxloc(T value) const;

    /**
     * @brief Find minimum value and its location
     * @tparam T Type of value (must be arithmetic type)
     * @param value Value from this process
     * @return MinLocResult containing min value and rank where it was found
     * @throws std::runtime_error if MPI operation fails
     * @throws std::invalid_argument if type T is not supported
     *
     * Finds the minimum value across all processes and returns both the
     * value and the rank of the process where it was found.
     */
    template<typename T>
    MinLocResult<T> minloc(T value) const;

    /**
     * @brief Broadcast value from root to all processes
     * @tparam T Type of value (must be arithmetic type)
     * @param value Value to broadcast (input on root, output on others)
     * @param root Rank of root process broadcasting the value
     * @throws std::runtime_error if MPI operation fails
     * @throws std::invalid_argument if type T is not supported
     *
     * Broadcasts a value from the root process to all other processes.
     * On the root process, value is the value to broadcast.
     * On other processes, value is overwritten with the broadcast value.
     */
    template<typename T>
    void broadcast(T& value, int root) const;

    /**
     * @brief Broadcast array from root to all processes
     * @tparam T Type of array elements (must be arithmetic type)
     * @param values Array to broadcast (input on root, output on others)
     * @param root Rank of root process broadcasting the array
     * @throws std::runtime_error if MPI operation fails
     * @throws std::invalid_argument if type T is not supported
     *
     * Broadcasts an array from the root process to all other processes.
     * On the root process, values is the array to broadcast.
     * On other processes, values is resized and overwritten with the broadcast array.
     */
    template<typename T>
    void broadcast_array(std::vector<T>& values, int root) const;

    /**
     * @brief Block until all processes reach this point
     * @throws std::runtime_error if barrier operation fails
     *
     * Calls MPI_Barrier to synchronize all processes in the communicator.
     * This is useful for coordinating timing or ensuring all processes
     * have completed a phase of computation.
     */
    void barrier() const;

    /**
     * @brief Get the communicator being used
     * @return MPI communicator
     */
    MPI_Comm communicator() const { return comm_; }

private:
    MPI_Comm comm_;  ///< MPI communicator for operations

    /**
     * @brief Get MPI_Op for a given ReductionOp
     * @param op Reduction operation
     * @return Corresponding MPI_Op
     * @throws std::invalid_argument if operation not supported
     */
    static MPI_Op get_mpi_op(ReductionOp op);

    /**
     * @brief Get MPI_Datatype for a given type
     * @tparam T Type to get datatype for
     * @return Corresponding MPI_Datatype
     * @throws std::invalid_argument if type not supported
     */
    template<typename T>
    static MPI_Datatype get_mpi_datatype();

    /**
     * @brief Get identity value for a reduction operation
     * @tparam T Type of value
     * @param op Reduction operation
     * @return Identity value for the operation
     */
    template<typename T>
    static T get_identity_value(ReductionOp op);
};

// Template method implementations

template<typename T>
T ReductionOps::reduce(T value, ReductionOp op, int root) const {
    MPI_Datatype datatype = get_mpi_datatype<T>();
    MPI_Op mpi_op = get_mpi_op(op);

    T result;
    int mpi_result = MPI_Reduce(&value, &result, 1, datatype, mpi_op, root, comm_);

    if (mpi_result != MPI_SUCCESS) {
        throw std::runtime_error("MPI_Reduce failed with error code: " + std::to_string(mpi_result));
    }

    return result;
}

template<typename T>
T ReductionOps::allreduce(T value, ReductionOp op) const {
    MPI_Datatype datatype = get_mpi_datatype<T>();
    MPI_Op mpi_op = get_mpi_op(op);

    T result;
    int mpi_result = MPI_Allreduce(&value, &result, 1, datatype, mpi_op, comm_);

    if (mpi_result != MPI_SUCCESS) {
        throw std::runtime_error("MPI_Allreduce failed with error code: " + std::to_string(mpi_result));
    }

    return result;
}

template<typename T>
T ReductionOps::scan(T value, ReductionOp op) const {
    MPI_Datatype datatype = get_mpi_datatype<T>();
    MPI_Op mpi_op = get_mpi_op(op);

    T result;
    int mpi_result = MPI_Scan(&value, &result, 1, datatype, mpi_op, comm_);

    if (mpi_result != MPI_SUCCESS) {
        throw std::runtime_error("MPI_Scan failed with error code: " + std::to_string(mpi_result));
    }

    return result;
}

template<typename T>
T ReductionOps::exscan(T value, ReductionOp op) const {
    MPI_Datatype datatype = get_mpi_datatype<T>();
    MPI_Op mpi_op = get_mpi_op(op);

    T result = get_identity_value<T>(op);
    int mpi_result = MPI_Exscan(&value, &result, 1, datatype, mpi_op, comm_);

    if (mpi_result != MPI_SUCCESS) {
        throw std::runtime_error("MPI_Exscan failed with error code: " + std::to_string(mpi_result));
    }

    return result;
}

template<typename T>
std::vector<T> ReductionOps::reduce_array(const std::vector<T>& values,
                                           ReductionOp op, int root) const {
    if (values.empty()) {
        return std::vector<T>();
    }

    MPI_Datatype datatype = get_mpi_datatype<T>();
    MPI_Op mpi_op = get_mpi_op(op);
    int count = static_cast<int>(values.size());

    std::vector<T> result(values.size());

    int mpi_result = MPI_Reduce(const_cast<T*>(values.data()), result.data(),
                                count, datatype, mpi_op, root, comm_);

    if (mpi_result != MPI_SUCCESS) {
        throw std::runtime_error("MPI_Reduce failed with error code: " + std::to_string(mpi_result));
    }

    // On non-root processes, return empty vector
    int rank;
    MPI_Comm_rank(comm_, &rank);
    if (rank != root) {
        return std::vector<T>();
    }

    return result;
}

template<typename T>
std::vector<T> ReductionOps::allreduce_array(const std::vector<T>& values,
                                              ReductionOp op) const {
    if (values.empty()) {
        return std::vector<T>();
    }

    MPI_Datatype datatype = get_mpi_datatype<T>();
    MPI_Op mpi_op = get_mpi_op(op);
    int count = static_cast<int>(values.size());

    std::vector<T> result(values.size());

    int mpi_result = MPI_Allreduce(const_cast<T*>(values.data()), result.data(),
                                   count, datatype, mpi_op, comm_);

    if (mpi_result != MPI_SUCCESS) {
        throw std::runtime_error("MPI_Allreduce failed with error code: " + std::to_string(mpi_result));
    }

    return result;
}

template<typename T>
MaxLocResult<T> ReductionOps::maxloc(T value) const {
    int rank;
    MPI_Comm_rank(comm_, &rank);

    struct {
        T value;
        int rank;
    } in_struct, out_struct;

    in_struct.value = value;
    in_struct.rank = rank;

    MPI_Datatype value_type = get_mpi_datatype<T>();

    // Create a struct datatype for {value, rank}
    int block_lengths[2] = {1, 1};
    MPI_Datatype types[2] = {value_type, MPI_INT};
    MPI_Aint displacements[2];

    MPI_Get_address(&in_struct.value, &displacements[0]);
    MPI_Get_address(&in_struct.rank, &displacements[1]);
    displacements[1] -= displacements[0];
    displacements[0] = 0;

    MPI_Datatype struct_type;
    MPI_Type_create_struct(2, block_lengths, displacements, types, &struct_type);
    MPI_Type_commit(&struct_type);

    int mpi_result = MPI_Allreduce(&in_struct, &out_struct, 1, struct_type, MPI_MAXLOC, comm_);

    MPI_Type_free(&struct_type);

    if (mpi_result != MPI_SUCCESS) {
        throw std::runtime_error("MPI_Allreduce failed with error code: " + std::to_string(mpi_result));
    }

    MaxLocResult<T> result;
    result.value = out_struct.value;
    result.rank = out_struct.rank;

    return result;
}

template<typename T>
MinLocResult<T> ReductionOps::minloc(T value) const {
    int rank;
    MPI_Comm_rank(comm_, &rank);

    struct {
        T value;
        int rank;
    } in_struct, out_struct;

    in_struct.value = value;
    in_struct.rank = rank;

    MPI_Datatype value_type = get_mpi_datatype<T>();

    // Create a struct datatype for {value, rank}
    int block_lengths[2] = {1, 1};
    MPI_Datatype types[2] = {value_type, MPI_INT};
    MPI_Aint displacements[2];

    MPI_Get_address(&in_struct.value, &displacements[0]);
    MPI_Get_address(&in_struct.rank, &displacements[1]);
    displacements[1] -= displacements[0];
    displacements[0] = 0;

    MPI_Datatype struct_type;
    MPI_Type_create_struct(2, block_lengths, displacements, types, &struct_type);
    MPI_Type_commit(&struct_type);

    int mpi_result = MPI_Allreduce(&in_struct, &out_struct, 1, struct_type, MPI_MINLOC, comm_);

    MPI_Type_free(&struct_type);

    if (mpi_result != MPI_SUCCESS) {
        throw std::runtime_error("MPI_Allreduce failed with error code: " + std::to_string(mpi_result));
    }

    MinLocResult<T> result;
    result.value = out_struct.value;
    result.rank = out_struct.rank;

    return result;
}

template<typename T>
void ReductionOps::broadcast(T& value, int root) const {
    MPI_Datatype datatype = get_mpi_datatype<T>();

    int mpi_result = MPI_Bcast(&value, 1, datatype, root, comm_);

    if (mpi_result != MPI_SUCCESS) {
        throw std::runtime_error("MPI_Bcast failed with error code: " + std::to_string(mpi_result));
    }
}

template<typename T>
void ReductionOps::broadcast_array(std::vector<T>& values, int root) const {
    int rank;
    MPI_Comm_rank(comm_, &rank);

    int size = 0;
    if (rank == root) {
        size = static_cast<int>(values.size());
    }

    // Broadcast the size first
    broadcast(size, root);

    // Resize the vector on non-root processes
    if (rank != root) {
        values.resize(size);
    }

    if (size > 0) {
        MPI_Datatype datatype = get_mpi_datatype<T>();

        int mpi_result = MPI_Bcast(values.data(), size, datatype, root, comm_);

        if (mpi_result != MPI_SUCCESS) {
            throw std::runtime_error("MPI_Bcast failed with error code: " + std::to_string(mpi_result));
        }
    }
}

// Template implementations for get_mpi_datatype

template<>
inline MPI_Datatype ReductionOps::get_mpi_datatype<int>() {
    return MPI_INT;
}

template<>
inline MPI_Datatype ReductionOps::get_mpi_datatype<long>() {
    return MPI_LONG;
}

template<>
inline MPI_Datatype ReductionOps::get_mpi_datatype<long long>() {
    return MPI_LONG_LONG;
}

template<>
inline MPI_Datatype ReductionOps::get_mpi_datatype<unsigned>() {
    return MPI_UNSIGNED;
}

template<>
inline MPI_Datatype ReductionOps::get_mpi_datatype<unsigned long>() {
    return MPI_UNSIGNED_LONG;
}

template<>
inline MPI_Datatype ReductionOps::get_mpi_datatype<unsigned long long>() {
    return MPI_UNSIGNED_LONG_LONG;
}

template<>
inline MPI_Datatype ReductionOps::get_mpi_datatype<float>() {
    return MPI_FLOAT;
}

template<>
inline MPI_Datatype ReductionOps::get_mpi_datatype<double>() {
    return MPI_DOUBLE;
}

template<>
inline MPI_Datatype ReductionOps::get_mpi_datatype<long double>() {
    return MPI_LONG_DOUBLE;
}

// Template implementations for get_identity_value

template<typename T>
inline T ReductionOps::get_identity_value(ReductionOp op) {
    switch (op) {
        case ReductionOp::SUM:
        case ReductionOp::LAND:
        case ReductionOp::LOR:
        case ReductionOp::LXOR:
            return static_cast<T>(0);
        case ReductionOp::PROD:
        case ReductionOp::BAND:
        case ReductionOp::BOR:
        case ReductionOp::BXOR:
            return static_cast<T>(1);
        case ReductionOp::MAX:
            return std::numeric_limits<T>::lowest();
        case ReductionOp::MIN:
            return std::numeric_limits<T>::max();
        default:
            return static_cast<T>(0);
    }
}

#endif // REDUCTION_OPS_HPP
