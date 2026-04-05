#ifndef CARTESIAN_TOPOLOGY_HPP
#define CARTESIAN_TOPOLOGY_HPP

#include <mpi.h>
#include <vector>
#include <stdexcept>
#include <cmath>

/**
 * @enum Direction
 * @brief Direction in Cartesian topology
 *
 * Used to specify which dimension to operate on.
 */
enum class Direction { X = 0, Y = 1 };

/**
 * @enum Shift
 * @brief Shift direction in Cartesian topology
 *
 * Used to specify which direction to shift coordinates.
 */
enum class Shift { Forward = 1, Backward = -1 };

/**
 * @struct NeighborInfo
 * @brief Information about neighboring processes in Cartesian topology
 *
 * Stores the ranks of the four neighbors (north, south, east, west).
 * Boundary processes have MPI_PROC_NULL as their neighbor rank.
 */
struct NeighborInfo {
    int south;  ///< Neighbor in -Y direction
    int north;  ///< Neighbor in +Y direction
    int west;   ///< Neighbor in -X direction
    int east;   ///< Neighbor in +X direction

    /**
     * @brief Construct NeighborInfo with MPI_PROC_NULL for all neighbors
     */
    NeighborInfo() : south(MPI_PROC_NULL), north(MPI_PROC_NULL),
                     west(MPI_PROC_NULL), east(MPI_PROC_NULL) {}
};

/**
 * @class CartesianTopology
 * @brief Manages MPI Cartesian topology for domain decomposition
 *
 * This class creates and manages a 2D Cartesian topology for MPI processes.
 * It provides functionality for:
 * - Automatic dimension calculation (closest to square)
 * - Coordinate mapping between ranks and positions
 * - Neighbor identification for communication patterns
 * - Boundary detection
 *
 * Features:
 * - RAII pattern: automatic cleanup on destruction
 * - Move-only: prevents accidental copies
 * - Exception-safe: handles errors gracefully
 * - Non-periodic boundaries (by default)
 */
class CartesianTopology {
public:
    /**
     * @brief Construct Cartesian topology with automatic dimension calculation
     * @param comm MPI communicator to create Cartesian topology from
     *
     * Automatically calculates optimal dimensions (closest to square) based
     * on the number of processes in the communicator.
     *
     * @throws std::runtime_error if topology creation fails
     */
    explicit CartesianTopology(MPI_Comm comm = MPI_COMM_WORLD);

    /**
     * @brief Construct Cartesian topology with specified dimensions
     * @param comm MPI communicator to create Cartesian topology from
     * @param dims Dimensions [dim_x, dim_y]
     *
     * Uses user-specified dimensions. The product of dimensions must equal
     * the number of processes in the communicator.
     *
     * @throws std::invalid_argument if dims.size() != 2
     * @throws std::invalid_argument if product of dims != number of processes
     * @throws std::runtime_error if topology creation fails
     */
    CartesianTopology(MPI_Comm comm, const std::vector<int>& dims);

    /**
     * @brief Destructor - free Cartesian communicator
     *
     * If this object owns the Cartesian communicator, it will be freed
     * using MPI_Comm_free.
     */
    ~CartesianTopology();

    // Delete copy constructor and copy assignment
    CartesianTopology(const CartesianTopology&) = delete;
    CartesianTopology& operator=(const CartesianTopology&) = delete;

    /**
     * @brief Move constructor
     * @param other The CartesianTopology to move from
     *
     * Transfers ownership of the Cartesian communicator.
     * The moved-from object becomes invalid and should not be used.
     */
    CartesianTopology(CartesianTopology&& other) noexcept;

    /**
     * @brief Move assignment operator
     * @param other The CartesianTopology to move from
     * @return Reference to this object
     *
     * Transfers ownership of the Cartesian communicator.
     * If this object currently owns a communicator, it is freed first.
     * The moved-from object becomes invalid and should not be used.
     */
    CartesianTopology& operator=(CartesianTopology&& other) noexcept;

    /**
     * @brief Get the Cartesian communicator
     * @return MPI communicator with Cartesian topology
     */
    MPI_Comm communicator() const { return cart_comm_; }

    /**
     * @brief Get the rank of the current process
     * @return Process rank in the Cartesian communicator
     */
    int rank() const { return rank_; }

    /**
     * @brief Get the total number of processes
     * @return Total number of processes in the Cartesian communicator
     */
    int size() const { return size_; }

    /**
     * @brief Get the dimensions of the topology
     * @return Vector [dim_x, dim_y]
     */
    const std::vector<int>& dims() const { return dims_; }

    /**
     * @brief Get the coordinates of the current process
     * @return Vector [coord_x, coord_y]
     */
    const std::vector<int>& coords() const { return coords_; }

    /**
     * @brief Get neighbor information
     * @return NeighborInfo structure with neighbor ranks
     */
    const NeighborInfo& neighbors() const { return neighbors_; }

    /**
     * @brief Get X dimension (number of processes in X direction)
     * @return Number of processes in X direction
     */
    int dim_x() const { return dims_[0]; }

    /**
     * @brief Get Y dimension (number of processes in Y direction)
     * @return Number of processes in Y direction
     */
    int dim_y() const { return dims_[1]; }

    /**
     * @brief Get X coordinate of current process
     * @return X coordinate (0 to dim_x-1)
     */
    int coord_x() const { return coords_[0]; }

    /**
     * @brief Get Y coordinate of current process
     * @return Y coordinate (0 to dim_y-1)
     */
    int coord_y() const { return coords_[1]; }

    /**
     * @brief Check if current process is on a boundary
     * @param dir Direction to check
     * @param shift Shift direction (Forward = +1, Backward = -1)
     * @return true if on boundary, false otherwise
     *
     * Examples:
     * - is_on_boundary(Direction::X, Shift::Backward) checks if on west boundary
     * - is_on_boundary(Direction::X, Shift::Forward) checks if on east boundary
     * - is_on_boundary(Direction::Y, Shift::Backward) checks if on south boundary
     * - is_on_boundary(Direction::Y, Shift::Forward) checks if on north boundary
     */
    bool is_on_boundary(Direction dir, Shift shift) const;

    /**
     * @brief Check if current process has a valid neighbor
     * @param dir Direction to check
     * @param shift Shift direction
     * @return true if neighbor has valid rank (not MPI_PROC_NULL)
     *
     * This is the opposite of is_on_boundary.
     */
    bool has_neighbor(Direction dir, Shift shift) const;

private:
    MPI_Comm cart_comm_;      ///< Cartesian communicator
    MPI_Comm parent_comm_;    ///< Parent communicator (for cleanup)
    int rank_;                ///< Process rank
    int size_;                ///< Total number of processes
    std::vector<int> dims_;   ///< Dimensions [dim_x, dim_y]
    std::vector<int> coords_; ///< Coordinates [coord_x, coord_y]
    NeighborInfo neighbors_;  ///< Neighbor ranks
    bool owns_comm_;          ///< True if this object owns cart_comm_

    /**
     * @brief Automatically compute optimal dimensions (closest to square)
     * @param num_procs Number of processes
     * @return Vector [dim_x, dim_y]
     *
     * Finds the dimensions that are closest to a square arrangement.
     * For example:
     * - 4 processes -> [2, 2]
     * - 6 processes -> [2, 3]
     * - 8 processes -> [2, 4]
     * - 12 processes -> [3, 4]
     */
    static std::vector<int> compute_optimal_dimensions(int num_procs);

    /**
     * @brief Compute neighbor ranks using MPI_Cart_shift
     *
     * Computes the ranks of the four neighbors (north, south, east, west)
     * and stores them in the neighbors_ structure.
     */
    void compute_neighbors();

    /**
     * @brief Validate that the Cartesian communicator is still valid
     * @throws std::runtime_error if communicator is MPI_COMM_NULL
     */
    void validate_communicator() const;
};

#endif // CARTESIAN_TOPOLOGY_HPP
