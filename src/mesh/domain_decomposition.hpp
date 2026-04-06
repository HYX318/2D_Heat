#ifndef DOMAIN_DECOMPOSITION_HPP
#define DOMAIN_DECOMPOSITION_HPP

#include "../mpi/cartesian_topology.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>

/**
 * @enum DecompositionStrategy
 * @brief Strategy for domain decomposition
 *
 * Different approaches to divide the global domain among MPI processes.
 */
enum class DecompositionStrategy {
    Automatic,     ///< Automatically automatically compute optimal decomposition
    Cartesian,     ///< Regular Cartesian grid
    Manual         ///< User-specified decomposition
};

/**
 * @struct Subdomain
 * @brief Information about a local subdomain
 *
 * Contains the dimensions and position of a subdomain in the global grid.
 * Includes coordinate conversion methods for mapping between global and local indices.
 */
struct Subdomain {
    size_t start_i;      ///< Global starting index in X direction
    size_t start_j;      ///< Global starting index in Y direction
    size_t nx;           ///< Local domain size in X direction
    size_t ny;           ///< Local domain size in Y direction
    size_t global_nx;    ///< Global domain size in X direction
    size_t global_ny;    ///< Global domain size in Y direction

    /**
     * @brief Convert global X index to local X index
     * @param global_i Global X index
     * @return Local X index (0 to nx-1)
     * @throws std::out_of_range if global_i is not in this subdomain
     */
    size_t global_to_local_i(size_t global_i) const {
        if (global_i < start_i || global_i >= start_i + nx) {
            throw std::out_of_range("Global X index is not in this subdomain");
        }
        return global_i - start_i;
    }

    /**
     * @brief Convert global Y index to local Y index
     * @param global_j Global Y index
     * @return Local Y index (0 to ny-1)
     * @throws std::out_of_range if global_j is not in this subdomain
     */
    size_t global_to_local_j(size_t global_j) const {
        if (global_j < start_j || global_j >= start_j + ny) {
            throw std::out_of_range("Global Y index is not in this subdomain");
        }
        return global_j - start_j;
    }
};

/**
 * @struct LoadBalanceStats
 * @brief Statistics about load balancing
 *
 * Provides information about how evenly the domain is distributed across processes.
 */
struct LoadBalanceStats {
    double min_load;            ///< Minimum load across all subdomains
    double max_load;            ///< Maximum load across all subdomains
    double avg_load;            ///< Average load across all subdomains
    double imbalance_ratio;     ///< Load imbalance ratio: (max - min) / avg

    LoadBalanceStats() : min_load(0.0), max_load(0.0), avg_load(0.0), imbalance_ratio(0.0) {}
};

/**
 * @class DomainDecomposition
 * @brief Manages domain decomposition for parallel computing
 *
 * This class divides a global 2D domain into subdomains for multiple MPI processes.
 * It supports:
 * - Automatic optimal decomposition (closest to square)
 * - User-specified Cartesian decomposition
 * - Manual non-uniform decomposition
 * - Coordinate transformation between global and local indices
 * - Boundary detection
 * - Load balance analysis
 *
 * Features:
 * - RAII pattern: automatic cleanup on destruction
 * - Move-only: prevents accidental copies
 * - Exception-safe: handles errors gracefully
 * - Flexible decomposition strategies
 */
class DomainDecomposition {
public:
    /**
     * @brief Construct with automatic decomposition
     * @param global_nx Global domain size in X direction
     * @param global_ny Global domain size in Y direction
     * @param topology Cartesian topology for MPI processes
     *
     * Automatically calculates optimal dimensions (closest to square) based
     * on the number of processes in the topology. Creates a regular Cartesian grid.
     *
     * @throws std::runtime_error if topology is invalid
     */
    DomainDecomposition(size_t global_nx, size_t global_ny,
                      const CartesianTopology& topology);

    /**
     * @brief Construct with specified Cartesian dimensions
     * @param global_nx Global domain size in X direction
     * @param global_ny Global domain size in Y direction
     * @param dim_x Number of processes in X direction
     * @param dim_y Number of processes in Y direction
     * @param topology Cartesian topology for MPI processes
     *
     * Uses user-specified dimensions. The product dim_x * dim_y must equal
     * the number of processes in the topology. Creates a regular Cartesian grid.
     *
     * @throws std::invalid_argument if dim_x * dim_y != topology.size()
     * @throws std::runtime_error if topology is invalid
     */
    DomainDecomposition(size_t global_nx, size_t global_ny,
                      size_t dim_x, size_t dim_y,
                      const CartesianTopology& topology);

    /**
     * @brief Construct with manual non-uniform decomposition
     * @param subdomains Vector of subdomain specifications
     * @param topology Cartesian topology for MPI processes
     *
     * Allows non-uniform decomposition where subdomains can have different sizes.
     * The number of subdomains must match the number of processes.
     *
     * @throws std::invalid_argument if subdomains.size() != topology.size()
     * @throws std::invalid_argument if subdomains don't cover the entire domain
     */
    DomainDecomposition(std::vector<Subdomain> subdomains,
                      const CartesianTopology& topology);

    /**
     * @brief Destructor
     */
    ~DomainDecomposition() = default;

    // Delete copy constructor and copy assignment
    DomainDecomposition(const DomainDecomposition&) = delete;
    DomainDecomposition& operator=(const DomainDecomposition&) = delete;

    /**
     * @brief Move constructor
     * @param other The DomainDecomposition to move from
     *
     * Transfers ownership of all data members.
     * The moved-from object becomes invalid and should not be used.
     */
    DomainDecomposition(DomainDecomposition&& other) noexcept = default;

    /**
     * @brief Move assignment operator
     * @param other The DomainDecomposition to move from
     * @return Reference to this object
     *
     * Transfers ownership of all data members.
     * The moved-from object becomes invalid and should not be used.
     */
    DomainDecomposition& operator=(DomainDecomposition&& other) noexcept = default;

    // ========== Accessors ==========

    /**
     * @brief Get the subdomain for the current process
     * @return Const reference to the current process's subdomain
     */
    const Subdomain& my_subdomain() const { return subdomains_[my_rank_]; }

    /**
     * @brief Get the rank of the current process
     * @return Process rank in the Cartesian communicator
     */
    size_t my_rank() const { return my_rank_; }

    /**
     * @brief Get all subdomains
     * @return Const reference to vector of all subdomains
     */
    const std::vector<Subdomain>& subdomains() const { return subdomains_; }

    /**
     * @brief Get global domain size in X direction
     * @return Global nx
     */
    size_t global_nx() const { return global_nx_; }

    /**
     * @brief Get global domain size in Y direction
     * @return Global ny
     */
    size_t global_ny() const { return global_ny_; }

    /**
     * @brief Get Cartesian dimension in X direction
     * @return Number of processes in X direction
     */
    size_t dim_x() const { return dim_x_; }

    /**
     * @brief Get Cartesian dimension in Y direction
     * @return Number of processes in Y direction
     */
    size_t dim_y() const { return dim_y_; }

    // ========== Decomposition Query ==========

    /**
     * @brief Check if global coordinates are in current process's subdomain
     * @param global_i Global X index
     * @param global_j Global Y index
     * @return true if coordinates are in this subdomain, false otherwise
     */
    bool is_in_my_domain(size_t global_i, size_t global_j) const;

    /**
     * @brief Find the rank that owns the given global coordinates
     * @param global_i Global X index
     * @param global_j Global Y index
     * @return Rank of the process that owns this cell
     * @throws std::out_of_range if coordinates are out of bounds
     */
    int find_owner_rank(size_t global_i, size_t global_j) const;

    /**
     * @brief Check if current process is on a boundary
     * @param dir Direction to check (X or Y)
     * @return true if on boundary in that direction, false otherwise
     *
     * For Direction::X: checks if on west or east boundary
     * For Direction::Y: checks if on south or north boundary
     */
    bool is_on_boundary(Direction dir) const;

    /**
     * @brief Check if current process has a neighbor in the given direction
     * @param dir Direction to check
     * @return true if neighbor exists, false otherwise
     *
     * This is the opposite of is_on_boundary.
     */
    bool has_neighbor(Direction dir) const;

    // ========== Load Analysis ==========

    /**
     * @brief Compute load balance statistics
     * @return LoadBalanceStats structure with load information
     *
     * Calculates the minimum, maximum, and average load across all subdomains,
     * as well as the imbalance ratio.
     */
    LoadBalanceStats compute_load_balance() const;

    // ========== Utility Methods ==========

    /**
     * @brief Compute optimal decomposition dimensions
     * @param num_procs Number of processes
     * @param global_nx Global domain size in X direction
     * @param global_ny Global domain size in Y direction
     * @return Pair [dim_x, dim_y] that is closest to a square arrangement
     *
     * Finds dimensions that minimize the aspect ratio difference and
     * balance the load as evenly as possible.
     */
    static std::pair<size_t, size_t> compute_optimal_dimensions(
        size_t num_procs, size_t global_nx, size_t global_ny);

private:
    size_t global_nx_;                    ///< Global domain size in X
    size_t global_ny_;                    ///< Global domain size in Y
    size_t dim_x_;                        ///< Number of processes in X direction
    size_t dim_y_;                        ///< Number of processes in Y direction
    size_t my_rank_;                      ///< Current process rank
    std::vector<Subdomain> subdomains_;  ///< All subdomains
    const CartesianTopology* topology_;   ///< Pointer to Cartesian topology

    /**
     * @brief Create regular Cartesian decomposition
     *
     * Divides the global domain into a regular Cartesian grid based on
     * the dimensions dim_x_ and dim_y_.
     */
    void create_cartesian_decomposition();

    /**
     * @brief Validate that all subdomains cover the global domain without overlap
     * @throws std::invalid_argument if validation fails
     */
    void validate_manual_decomposition() const;

    /**
     * @brief Validate that Cartesian topology is still valid
     * @throws std::runtime_error if topology is invalid
     */
    void validate_topology() const;
};

#endif // DOMAIN_DECOMPOSITION_HPP
