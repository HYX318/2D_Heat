/**
 * @file mesh2d.hpp
 * @brief 2D rectangular mesh class with ghost cell support
 *
 * This class manages 2D rectangular grids with support for ghost cells,
 * boundary conditions, and MPI domain decomposition.
 */

#ifndef MESH2D_HPP
#define MESH2D_HPP

#include <cstddef>
#include <functional>
#include <memory>
#include <stdexcept>
#include <cmath>
#include "../utils/array2d.hpp"
#include "../mpi/cartesian_topology.hpp"
#include "../mpi/ghost_cell_exchange.hpp"
#include "../utils/logger.hpp"

/**
 * @enum Direction
 * @brief Cardinal directions for boundary conditions and ghost cell operations
 */
enum class Direction {
    North = 0,  ///< +Y direction
    South = 1,  ///< -Y direction
    East = 2,   ///< +X direction
    West = 3    ///< -X direction
};

/**
 * @class Mesh2D
 * @brief 2D rectangular mesh with ghost cell support
 *
 * This class manages a 2D rectangular grid with support for:
 * - Ghost cells for MPI domain decomposition
 * - Dirichlet boundary conditions
 * - Coordinate transformations (local/global, grid/physical)
 * - Numerical operations (Laplacian, norms, arithmetic)
 * - Efficient ghost cell exchange
 *
 * Grid layout with ghost cells:
 * ```
 * Row 0      : South ghost row
 * Row [1..ny]: Interior rows
 * Row ny+1   : North ghost row
 *
 * Col 0      : West ghost column
 * Col [1..nx]: Interior columns
 * Col nx+1   : East ghost column
 * ```
 *
 * Features:
 * - RAII pattern for automatic resource management
 * - Move semantics for efficient transfers
 * - Exception-safe operations
 * - MPI-aware ghost cell exchange
 * - Efficient 5-point stencil operations
 */
class Mesh2D {
public:
    /**
     * @brief Basic constructor - creates local mesh without ghost cells
     * @param nx Number of interior points in x-direction
     * @param ny Number of interior points in y-direction
     * @param lx Physical length in x-direction
     * @param ly Physical length in y-direction
     * @throws std::invalid_argument if nx or ny is 0
     */
    Mesh2D(size_t nx, size_t ny, double lx = 1.0, double ly = 1.0);

    /**
     * @brief Constructor with MPI - creates mesh with ghost cells
     * @param nx Global number of interior points in x-direction
     * @param ny Global number of interior points in y-direction
     * @param lx Physical length in x-direction
     * @param ly Physical length in y-direction
     * @param topology Cartesian topology for domain decomposition
     * @throws std::invalid_argument if nx or ny is 0
     * @throws std::runtime_error if memory allocation fails
     */
    Mesh2D(size_t nx, size_t ny, double lx, double ly,
           const CartesianTopology& topology);

    /**
     * @brief Copy constructor - performs deep copy
     * @param other Mesh to copy from
     */
    Mesh2D(const Mesh2D& other);

    /**
     * @brief Move constructor - efficient transfer of ownership
     * @param other Mesh to move from
     */
    Mesh2D(Mesh2D&& other) noexcept;

    /**
     * @brief Destructor - automatically handles cleanup
     */
    ~Mesh2D();

    /**
     * @brief Copy assignment operator
     * @param other Mesh to copy from
     * @return Reference to this mesh
     */
    Mesh2D& operator=(const Mesh2D& other);

    /**
     * @brief Move assignment operator
     * @param other Mesh to move from
     * @return Reference to this mesh
     */
    Mesh2D& operator=(Mesh2D&& other) noexcept;

    /**
     * @brief Element access (mutable) - accounts for ghost cell offset
     * @param i Row index (0 to ny-1 for interior)
     * @param j Column index (0 to nx-1 for interior)
     * @return Reference to element at (i, j)
     * @throws std::out_of_range if indices are out of bounds
     */
    double& operator()(size_t i, size_t j);

    /**
     * @brief Element access (const) - accounts for ghost cell offset
     * @param i Row index (0 to ny-1 for interior)
     * @param j Column index (0 to nx-1 for interior)
     * @return Const reference to element at (i, j)
     * @throws std::out_of_range if indices are out of bounds
     */
    const double& operator()(size_t i, size_t j) const;

    /**
     * @brief Direct element access to underlying array (mutable)
     * @param i Array row index (includes ghost cells)
     * @param j Array column index (includes ghost cells)
     * @return Reference to element at (i, j)
     * @throws std::out_of_range if indices are out of bounds
     */
    double& at(size_t i, size_t j);

    /**
     * @brief Direct element access to underlying array (const)
     * @param i Array row index (includes ghost cells)
     * @param j Array column index (includes ghost cells)
     * @return Const reference to element at (i, j)
     * @throws std::out_of_range if indices are out of bounds
     */
    const double& at(size_t i, size_t j) const;

    /**
     * @brief Get number of interior points in x-direction
     * @return nx
     */
    size_t nx() const noexcept { return nx_; }

    /**
     * @brief Get number of interior points in y-direction
     * @return ny
     */
    size_t ny() const noexcept { return ny_; }

    /**
     * @brief Get total number of points in x-direction (including ghost cells)
     * @return nx + 2 (if ghost cells enabled) or nx
     */
    size_t total_nx() const noexcept { return total_nx_; }

    /**
     * @brief Get total number of points in y-direction (including ghost cells)
     * @return ny + 2 (if ghost cells enabled) or ny
     */
    size_t total_ny() const noexcept { return total_ny_; }

    /**
     * @brief Get physical length in x-direction
     * @return lx
     */
    double lx() const noexcept { return lx_; }

    /**
     * @brief Get physical length in y-direction
     * @return ly
     */
    double ly() const noexcept { return ly_; }

    /**
     * @brief Get grid spacing in x-direction
     * @return hx = lx / (nx - 1) or lx / nx depending on interpretation
     */
    double hx() const noexcept { return hx_; }

    /**
     * @brief Get grid spacing in y-direction
     * @return hy = ly / (ny - 1) or ly / ny depending on interpretation
     */
    double hy() const noexcept { return hy_; }

    /**
     * @brief Check if mesh has ghost cells
     * @return true if ghost cells are enabled
     */
    bool has_ghost_cells() const noexcept { return has_ghost_cells_; }

    /**
     * @brief Convert global x-index to local x-index
     * @param global_i Global index in x-direction
     * @return Local index
     * @throws std::out_of_range if global index is out of local range
     */
    size_t global_to_local_x(size_t global_i) const;

    /**
     * @brief Convert global y-index to local y-index
     * @param global_j Global index in y-direction
     * @return Local index
     * @throws std::out_of_range if global index is out of local range
     */
    size_t global_to_local_y(size_t global_j) const;

    /**
     * @brief Convert local x-index to global x-index
     * @param local_i Local index in x-direction
     * @return Global index
     */
    size_t local_to_global_x(size_t local_i) const;

    /**
     * @brief Convert local y-index to global y-index
     * @param local_j Local index in y-direction
     * @return Global index
     */
    size_t local_to_global_y(size_t local_j) const;

    /**
     * @brief Get physical x-coordinate at grid point
     * @param i Grid index in x-direction
     * @return Physical x-coordinate
     * @throws std::out_of_range if index is out of bounds
     */
    double x_coord(size_t i) const;

    /**
     * @brief Get physical y-coordinate at grid point
     * @param j Grid index in y-direction
     * @return Physical y-coordinate
     * @throws std::out_of_range if index is out of bounds
     */
    double y_coord(size_t j) const;

    /**
     * @brief Get physical coordinates at grid point
     * @param i Grid index in x-direction
     * @param j Grid index in y-direction
     * @return Pair (x, y) of physical coordinates
     */
    std::pair<double, double> coord(size_t i, size_t j) const;

    /**
     * @brief Apply constant Dirichlet boundary condition
     * @param value Boundary value
     * @param t Time parameter (for time-dependent BCs)
     */
    void apply_dirichlet_bc(double value = 0.0, double t = 0.0);

    /**
     * @brief Apply function-based boundary condition
     * @param bc_func Boundary condition function f(x, y, t)
     * @param t Time parameter
     */
    void apply_bc(std::function<double(double x, double y, double t)> bc_func, double t = 0.0);

    /**
     * @brief Apply constant boundary condition at specific direction
     * @param dir Direction to apply BC to
     * @param value Boundary value
     */
    void apply_bc_at_direction(Direction dir, double value);

    /**
     * @brief Compute Laplacian using 5-point stencil
     * @param result Output array for Laplacian
     * @throws std::invalid_argument if result size doesn't match
     *
     * Uses central difference:
     * ∇²u = (u_{i+1,j} - 2u_{i,j} + u_{i-1,j})/hx² + (u_{i,j+1} - 2u_{i,j} + u_{i,j-1})/hy²
     */
    void compute_laplacian(utils::Array2D& result) const;

    /**
     * @brief Fill all interior points with specified value
     * @param value Fill value
     */
    void fill(double value);

    /**
     * @brief Copy data from another mesh
     * @param other Source mesh
     * @throws std::invalid_argument if dimensions don't match
     */
    void copy_from(const Mesh2D& other);

    /**
     * @brief Scale all interior points by factor
     * @param factor Scaling factor
     */
    void scale(double factor);

    /**
     * @brief Add another mesh to this mesh
     * @param other Mesh-mesh to add
     * @throws std::invalid_argument if dimensions don't match
     */
    void add(const Mesh2D& other);

    /**
     * @brief Subtract another mesh from this mesh
     * @param other Mesh to subtract
     * @throws std::invalid_argument if dimensions don't match
     */
    void subtract(const Mesh2D& other);

    /**
     * @brief Compute L2 norm (Euclidean norm)
     * @return sqrt(sum of squares of all interior elements)
     */
    double l2_norm() const;

    /**
     * @brief Compute L-infinity norm (maximum absolute value)
     * @return Maximum absolute value in interior elements
     */
    double linfty_norm() const;

    /**
     * @brief Get maximum value
     * @return Maximum value in interior elements
     */
    double max() const;

    /**
     * @brief Get minimum value
     * @return Minimum value in interior elements
     */
    double min() const;

    /**
     * @brief Exchange ghost cells with neighboring processes
     * @throws std::runtime_error if MPI communication fails or ghost cells not enabled
     */
    void exchange_ghost_cells();

    /**
     * @brief Set all boundary values (including ghost cells)
     * @param value Boundary value
     */
    void set_boundary_values(double value);

    /**
     * @brief Get reference to underlying data array
     * @return Reference to Array2D
     */
    utils::Array2D& data() noexcept { return data_; }

    /**
     * @brief Get const reference to underlying data array
     * @return Const reference to Array2D
     */
    const utils::Array2D& data() const noexcept { return data_; }

    /**
     * @brief Get pointer to Cartesian topology (may be nullptr)
     * @return Pointer to topology or nullptr if not using MPI
     */
    const CartesianTopology* topology() const noexcept { return topology_; }

    /**
     * @brief Get pointer to ghost cell exchange (may be nullptr)
     * @return Pointer to ghost cell exchange or nullptr if not using MPI
     */
    GhostCellExchange* ghost_exchange() const noexcept { return ghost_exchange_; }

private:
    // Grid dimensions
    size_t nx_;              ///< Interior points in x
    size_t ny_;              ///< Interior points in y
    size_t total_nx_;        ///< Total points in x (including ghost)
    size_t total_ny_;        ///< Total points in y (including ghost)

    // Physical dimensions
    double lx_;              ///< Physical length in x
    double ly_;              ///< Physical length in y
    double hx_;              ///< Grid spacing in x
    double hy_;              ///< Grid spacing in y

    // Ghost cell support
    bool has_ghost_cells_;  ///< True if ghost cells enabled

    // Data storage
    utils::Array2D data_;    ///< Array2D for data storage

    // MPI support
    const CartesianTopology* topology_;           ///< Pointer to topology (not owned)
    std::unique_ptr<GhostCellExchange> ghost_exchange_;  ///< Ghost cell exchange

    // Offset for ghost cells
    size_t offset_x_;        ///< Offset for x-index (1 if ghost cells, 0 otherwise)
    size_t offset_y_;        ///< Offset for y-index (1 if ghost cells, 0 otherwise)

    /**
     * @brief Validate indices are within interior bounds
     * @param i Row index
     * @param j Column index
     * @throws std::out_of_range if indices are invalid
     */
    void check_interior_bounds(size_t i, size_t j) const;

    /**
     * @brief Validate indices are within total bounds (including ghost)
     * @param i Row index
     * @param j Column index
     * @throws std::out_of_range if indices are invalid
     */
    void check_total_bounds(size_t i, size_t j) const;

    /**
     * @brief Validate index is within x bounds
     * @param i Index in x-direction
     * @throws std::out_of_range if index is invalid
     */
    void check_x_bounds(size_t i) const;

    /**
     * @brief Validate index is within y bounds
     * @param j Index in y-direction
     * @throws std::out_of_range if index is invalid
     */
    void check_y_bounds(size_t j) const;

    /**
     * @brief Compute grid spacings
     */
    void compute_spacings();

    /**
     * @brief Initialize ghost cell exchange
     * @param topology Cartesian topology reference
     */
    void initialize_ghost_exchange(const CartesianTopology& topology);
};

#endif // MESH2D_HPP
