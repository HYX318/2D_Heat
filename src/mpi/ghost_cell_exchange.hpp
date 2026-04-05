#ifndef GHOST_CELL_EXCHANGE_HPP
#define GHOST_CELL_EXCHANGE_HPP

#include <mpi.h>
#include "../utils/array2d.hpp"
#include "cartesian_topology.hpp"
#include <stdexcept>
#include <cassert>

/**
 * @namespace ghost
 * @brief Direction constants for ghost cell exchange
 */
namespace ghost {
    constexpr int SOUTH = 0;
    constexpr int NORTH = 1;
    constexpr int WEST = 2;
    constexpr int EAST = 3;
}

/**
 * @class GhostCellExchange
 * @brief Manages ghost cell exchange for MPI domain decomposition
 *
 * This class handles the exchange of ghost cells between neighboring
 * processes in a 2D Cartesian topology. It provides both synchronous
 * and asynchronous exchange operations using custom MPI data types.
 *
 * Features:
 * - Custom MPI data types for efficient row/column communication
 * - Synchronous exchange using MPI_Sendrecv (deadlock-safe)
 * - Asynchronous exchange using non-blocking communication
 * - Exception-safe RAII pattern
 * - Automatic boundary process detection
 */
class GhostCellExchange {
public:
    /**
     * @brief Construct GhostCellExchange with domain dimensions and topology
     * @param nx Number of interior columns (excluding ghost cells)
     * @param ny Number of interior rows (excluding ghost cells)
     * @param topology Cartesian topology for neighbor information
     * @throws std::invalid_argument if dimensions are invalid
     * @throws std::runtime_error if MPI data type creation fails
     *
     * Creates custom MPI data types for efficient communication of
     * rows and columns. The array layout with ghost cells is:
     *
     * Rows:    [0] = South ghost, [1..ny] = Interior, [ny+1] = North ghost
     * Columns: [0] = West ghost,  [1..nx] = Interior, [nx+1] = East ghost
     */
    GhostCellExchange(int nx, int ny, const CartesianTopology& topology);

    /**
     * @brief Destructor - free custom MPI data types
     *
     * Exception-safe: catches and logs any errors during cleanup.
     */
    ~GhostCellExchange();

    // Delete copy constructor and copy assignment
    GhostCellExchange(const GhostCellExchange&) = delete;
    GhostCellExchange& operator=(const GhostCellExchange&) = delete;

    /**
     * @brief Move constructor
     * @param other The GhostCellExchange to move from
     *
     * Transfers ownership of MPI data types.
     * The moved-from object becomes invalid.
     */
    GhostCellExchange(GhostCellExchange&& other) noexcept;

    /**
     * @brief Move assignment operator
     * @param other The GhostCellExchange to move from
     * @return Reference to this object
     *
     * Transfers ownership of MPI data types.
     * The moved-from object becomes invalid.
     */
    GhostCellExchange& operator=(GhostCellExchange&& other) noexcept;

    /**
     * @brief Perform synchronous ghost cell exchange
     * @param array Array2D to exchange ghost cells for
     * @throws std::invalid_argument if array size doesn't match
     * @throws std::runtime_error if MPI communication fails
     *
     * Exchange pattern (using MPI_Sendrecv to avoid deadlock):
     * 1. Send bottom row to South, receive   from North
     * 2. Send top row    to North, receive   from South
     * 3. Send left column to West,  receive   from East
     * 4. Send right column to East,  receive   from West
     *
     * Layout:
     * ```
     * ┌─────────────────────────────────┐
     * │  North ghost row (receive)     │
     * │  ┌───────────────────────────┐ │
     * │  │    Interior Domain       │ │
     * │  └───────────────────────────┘ │
     * │  South ghost row (receive)     │
     * └─────────────────────────────────┘
     * ```
     */
    void exchange(utils::Array2D& array);

    /**
     * @brief Start asynchronous ghost cell exchange
     * @param array Array2D to exchange ghost cells for
     * @param requests Array of 8 MPI_Request objects for communication
     * @throws std::invalid_argument if array size doesn't match
     * @throws std::runtime_error if MPI communication fails
     *
     * Initiates non-blocking communication for all ghost cell exchanges.
     * The requests array must have space for 8 MPI_Request objects:
     *
     * - requests[0]: Irecv from North  (for south boundary)
     * - requests[1]: Isend to   South  (from south boundary)
     * - requests[2]: Irecv from South  (for north boundary)
     * - requests[3]: Isend to   North  (from north boundary)
     * - requests[4]: Irecv from East   (for west boundary)
     * - requests[5]: Isend to   West   (from west boundary)
     * - requests[6]: Irecv from West   (for east boundary)
     * - requests[7]: Isend to   East   (from east boundary)
     *
     * After calling this function, call wait_all() to wait for completion.
     */
    void exchange_async(utils::Array2D& array, MPI_Request* requests);

    /**
     * @brief Get custom MPI data type for rows
     * @return MPI_Datatype for row communication
     */
    MPI_Datatype row_type() const { return row_type_; }

    /**
     * @brief Get custom MPI data type for columns
     * @return MPI_Datatype for column communication
     */
    MPI_Datatype col_type() const { return col_type_; }

    /**
     * @brief Validate array size matches expected dimensions
     * @param array Array to validate
     * @return true if size is correct, false otherwise
     */
    bool validate_array_size(const utils::Array2D& array) const;

    /**
     * @brief Wait for all asynchronous exchange requests to complete
     * @param requests Array of 8 MPI_Request objects from exchange_async
     * @throws std::runtime_error if MPI_Waitall fails
     *
     * Waits for all 8 requests to complete and frees the request objects.
     * The requests array will contain MPI_REQUEST_NULL after completion.
     */
    void wait_all(MPI_Request* requests) const;

    /**
     * @brief Get number of interior columns (excluding ghost cells)
     * @return Interior column count
     */
    int nx() const { return nx_; }

    /**
     * @brief Get number of interior rows (excluding ghost cells)
     * @return Interior row count
     */
    int ny() const { return ny_; }

    /**
     * @brief Get reference to Cartesian topology
     * @return Const reference to topology
     */
    const CartesianTopology& topology() const { return topology; }

    /**
     * @brief Get neighbor information
     * @return Const reference to neighbor info
     */
    const NeighborInfo& neighbors() const { return topology.neighbors(); }

private:
    int nx_;                           ///< Interior columns
    int ny_;                           ///< Interior rows
    const CartesianTopology& topology; ///< Reference to topology
    MPI_Datatype row_type_;            ///< Custom type for rows
    MPI_Datatype col_type_;            ///< Custom type for columns
    bool owns_types_;                  ///< True if we own the MPI types

    /**
     * @brief Create custom MPI data types
     * @throws std::runtime_error if type creation fails
     *
     * Creates:
     * - row_type_: contiguous type for a row (nx+2 elements)
     * - col_type_: vector type for a column (ny+2 elements, stride nx+2)
     */
    void create_datatypes();

    /**
     * @brief Free custom MPI data types
     *
     * Exception-safe: catches and logs any errors.
     */
    void free_datatypes();

    /**
     * @brief Validate dimensions are valid
     * @throws std::invalid_argument if dimensions are invalid
     */
    void validate_dimensions() const;
};

#endif // GHOST_CELL_EXCHANGE_HPP
