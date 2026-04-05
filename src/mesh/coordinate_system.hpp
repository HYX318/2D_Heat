/**
 * @file coordinate_system.hpp
 * @brief Coordinate system, grid generation, and coordinate transformations
 *
 * This file provides a comprehensive CoordinateSystem class for handling
 * various grid types including uniform, non-uniform, and stretched grids
 * commonly used in computational fluid dynamics and heat equation solvers.
 */

#ifndef COORDINATE_SYSTEM_HPP
#define COORDINATE_SYSTEM_HPP

#include <cstddef>
#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include "../utils/array2d.hpp"

namespace mesh {

/**
 * @enum GridType
 * @brief Different types of grid distributions
 */
enum class GridType {
    Uniform,        ///< Uniform grid spacing
    NonUniform,     ///< Non-uniform grid spacing
    Stretched       ///< Stretched grid (for boundary layers)
};

/**
 * @enum CoordinateSystemType
 * @brief Different coordinate systems
 */
enum class CoordinateSystemType {
    Cartesian,       ///< Cartesian coordinates (x, y)
    Polar,          ///< Polar coordinates (r, θ)
    Spherical        ///< Spherical coordinates (r, θ, φ)
};

/**
 * @struct GridParams
 * @brief Parameters for grid generation
 */
struct GridParams {
    size_t nx;                                             ///< X direction number of points
    size_t ny;                                             ///< Y direction number of points
    double x_min;                                          ///< X minimum value
    double x_max;                                          ///< X maximum value
    double y_min;                                          ///< Y minimum value
    double y_max;                                          ///< Y maximum value
    GridType type;                                         ///< Grid type

    // Non-uniform grid parameters (optional)
    std::function<double(double)> x_distribution;          ///< X coordinate distribution function
    std::function<double(double)> y_distribution;          ///< Y coordinate distribution function

    /**
     * @brief Default constructor with sensible defaults
     */
    GridParams()
        : nx(100), ny(100),
          x_min(0.0), x_max(1.0),
          y_min(0.0), y_max(1.0),
          type(GridType::Uniform),
          x_distribution(nullptr),
          y_distribution(nullptr) {}
};

/**
 * @class CoordinateSystem
 * @brief Handles coordinate system, grid generation, and coordinate transformations
 *
 * This class provides:
 * - Uniform and non-uniform grid generation
 * - Coordinate indexing and lookup
 * - Coordinate transformations
 * - Grid quality metrics
 * - Efficient coordinate calculations
 *
 * Example usage:
 * @code
 * // Uniform grid
 * CoordinateSystem coords(100, 100, 0.0, 1.0, 0.0, 1.0);
 *
 * // Access coordinates
 * double x = coords.x(50);    // Center x coordinate
 * double y = coords.y(50);    // Center y coordinate
 *
 * // Non-uniform stretched grid
 * CoordinateSystem coords_stretched(100, 100);
 * coords_stretched.set_x_distribution([](double xi) {
 *     return CoordinateSystem::stretching_function(xi, 1.2);
 * });
 *
 * // Generate coordinate mesh
 * utils::Array2D x_mesh(100, 100);
 * utils::Array2D y_mesh(100, 100);
 * coords.generate_coordinate_mesh(x_mesh, y_mesh);
 * @endcode
 */
class CoordinateSystem {
public:
    /**
     * @brief Constructor for uniform grid
     * @param nx Number of points in X direction
     * @param ny Number of points in Y direction
     * @param x_min Minimum X value
     * @param x_max Maximum X value
     * @param y_min Minimum Y value
     * @param y_max Maximum Y value
     * @throws std::invalid_argument if nx or ny is 0, or if ranges are invalid
     */
    CoordinateSystem(size_t nx, size_t ny,
                    double x_min = 0.0, double x_max = 1.0,
                    double y_min = 0.0, double y_max = 1.0);

    /**
     * @brief Constructor from GridParams
     * @param params Grid parameters
     * @throws std::invalid_argument if parameters are invalid
     */
    explicit CoordinateSystem(const GridParams& params);

    /**
     * @brief Copy constructor
     * @param other CoordinateSystem to copy from
     */
    CoordinateSystem(const CoordinateSystem& other) = default;

    /**
     * @brief Copy assignment operator
     * @param other CoordinateSystem to copy from
     * @return Reference to this object
     */
    CoordinateSystem& operator=(const CoordinateSystem& other) = default;

    /**
     * @brief Move constructor
     * @param other CoordinateSystem to move from
     */
    CoordinateSystem(CoordinateSystem&& other) noexcept = default;

    /**
     * @brief Move assignment operator
     * @param other CoordinateSystem to move from
     * @return Reference to this object
     */
    CoordinateSystem& operator=(CoordinateSystem&& other) noexcept = default;

    /**
     * @brief Destructor
     */
    ~CoordinateSystem() = default;

    // ========== Grid Information ==========

    /**
     * @brief Get number of points in X direction
     * @return nx
     */
    size_t nx() const noexcept { return nx_; }

    /**
     * @brief Get number of points in Y direction
     * @return ny
     */
    size_t ny() const noexcept { return ny_; }

    /**
     * @brief Get minimum X value
     * @return x_min
     */
    double x_min() const noexcept { return x_min_; }

    /**
     * @brief Get maximum X value
     * @return x_max
     */
    double x_max() const noexcept { return x_max_; }

    /**
     * @brief Get minimum Y value
     * @return y_min
     */
    double y_min() const noexcept { return y_min_; }

    /**
     * @brief Get maximum Y value
     * @return y_max
     */
    double y_max() const noexcept { return y_max_; }

    /**
     * @brief Get domain length in X direction
     * @return x_max - x_min
     */
    double lx() const noexcept { return x_max_ - x_min_; }

    /**
     * @brief Get domain length in Y direction
     * @return y_max - y_min
     */
    double ly() const noexcept { return y_max_ - y_min_; }

    /**
     * @brief Get average grid spacing in X direction
     * @return (x_max - x_min) / (nx - 1)
     */
    double hx() const noexcept { return lx() / (nx_ - 1); }

    /**
     * @brief Get average grid spacing in Y direction
     * @return (y_max - y_min) / (ny - 1)
     */
    double hy() const noexcept { return ly() / (ny_ - 1); }

    /**
     * @brief Get grid type
     * @return Grid type
     */
    GridType grid_type() const noexcept { return grid_type_; }

    // ========== Coordinate Calculation ==========

    /**
     * @brief Get X coordinate from index
     * @param i X index (0-based)
     * @return X coordinate
     * @throws std::out_of_range if i >= nx
     */
    double x(size_t i) const;

    /**
     * @brief Get Y coordinate from index
     * @param j Y index (0-based)
     * @return Y coordinate
     * @throws std::out_of_range if j >= ny
     */
    double y(size_t j) const;

    /**
     * @brief Get nearest index from X coordinate
     * @param x X coordinate
     * @return Nearest index (clamped to [0, nx-1])
     */
    size_t i_at(double x) const;

    /**
     * @brief Get nearest index from Y coordinate
     * @param y Y coordinate
     * @return Nearest index (clamped to [0, ny-1])
     */
    size_t j_at(double y) const;

    /**
     * @brief Get coordinate pair from indices
     * @param i X index
     * @param j Y index
     * @return Pair (x, y)
     */
    std::pair<double, double> coord(size_t i, size_t j) const;

    /**
     * @brief Get grid spacing at index i in X direction
     * @param i X index
     * @return Grid spacing (may vary for non-uniform grids)
     */
    double hx_at(size_t i) const;

    /**
     * @brief Get grid spacing at index j in Y direction
     * @param j Y index
     * @return Grid spacing (may vary for non-uniform grids)
     */
    double hy_at(size_t j) const;

    // ========== Coordinate Transformations ==========

    /**
     * @brief Get normalized X coordinate from [0, 1]
     * @param i X index
     * @return Normalized coordinate in [0, 1]
     */
    double xi(size_t i) const;

    /**
     * @brief Get normalized Y coordinate from [0, 1]
     * @param j Y index
     * @return Normalized coordinate in [0, 1]
     */
    double eta(size_t j) const;

    /**
     * @brief Convert physical X coordinate to normalized coordinate
     * @param x Physical X coordinate
     * @return Normalized coordinate in [0, 1]
     */
    double xi_from_x(double x) const;

    /**
     * @brief Convert physical Y coordinate to normalized coordinate
     * @param y Physical Y coordinate
     * @return Normalized coordinate in [0, 1]
     */
    double eta_from_y(double y) const;

    /**
     * @brief Convert normalized coordinate to physical X coordinate
     * @param xi Normalized coordinate in [0, 1]
     * @return Physical X coordinate
     */
    double x_from_xi(double xi) const;

    /**
     * @brief Convert normalized coordinate to physical Y coordinate
     * @param eta Normalized coordinate in [0, 1]
     * @return Physical Y coordinate
     */
    double y_from_eta(double eta) const;

    // ========== Grid Generation ==========

    /**
     * @brief Generate X coordinates array
     * @return Vector of X coordinates
     */
    std::vector<double> generate_x_coords() const;

    /**
     * @brief Generate Y coordinates array
     * @return Vector of Y coordinates
     */
    std::vector<double> generate_y_coords() const;

    /**
     * @brief Generate coordinate mesh
     * @param x_mesh Output X coordinate mesh (must match dimensions)
     * @param y_mesh Output Y coordinate mesh (must match dimensions)
     * @throws std::invalid_argument if dimensions don't match
     */
    void generate_coordinate_mesh(utils::Array2D& x_mesh,
                                  utils::Array2D& y_mesh) const;

    // ========== Non-Uniform Grid Support ==========

    /**
     * @brief Set non-uniform X distribution function
     * @param dist Function mapping normalized [0,1] to physical coordinate
     */
    void set_x_distribution(std::function<double(double)> dist);

    /**
     * @brief Set non-uniform Y distribution function
     * @param dist Function mapping normalized [0,1] to physical coordinate
     */
    void set_y_distribution(std::function<double(double)> dist);

    /**
     * @brief Stretching function for boundary layer grids
     * @param xi Normalized coordinate in [0, 1]
     * @param beta Stretching parameter (beta > 1 clusters points near 0)
     * @return Stretched coordinate in [0, 1]
     *
     * This function implements the standard stretching function:
     * ξ_stretched = (1 - β^ξ) / (1 - β)
     *
     * For beta > 1, points are clustered near ξ = 0
     * For beta < 1, points are clustered near ξ = 1
     */
    static double stretching_function(double xi, double beta);

    // ========== Validation and Utilities ==========

    /**
     * @brief Check if coordinate is inside domain
     * @param x X coordinate
     * @param y Y coordinate
     * @return true if inside domain (inclusive)
     */
    bool is_inside(double x, double y) const;

    /**
     * @brief Check if coordinate is inside domain with tolerance
     * @param x X coordinate
     * @param y Y coordinate
     * @param tolerance Tolerance for boundary checking
     * @return true if inside domain (with tolerance)
     */
    bool is_inside(double x, double y, double tolerance) const;

    /**
     * @brief Compute grid quality metric
     * @return Grid quality (closest to 1 is optimal)
     *
     * For uniform grids, quality = 1.0
     * For non-uniform grids, quality measures spacing variation
     */
    double grid_quality() const;

    /**
     * @brief Compute aspect ratio of grid
     * @return Aspect ratio (lx / ly)
     */
    double aspect_ratio() const;

private:
    size_t nx_;                                            ///< Number of points in X
    size_t ny_;                                            ///< Number of points in Y
    double x_min_;                                         ///< Minimum X value
    double x_max_;                                         ///< Maximum X value
    double y_min_;                                         ///< Minimum Y value
    double y_max_;                                         ///< Maximum Y value
    GridType grid_type_;                                   ///< Grid type

    // Distribution functions for non-uniform grids
    std::function<double(double)> x_distribution_;        ///< X distribution function
    std::function<double(double)> y_distribution_;        ///< Y distribution function

    // Cached coordinate arrays for non-uniform grids (lazily computed)
    mutable std::vector<double> x_coords_cache_;          ///< Cached X coordinates
    mutable std::vector<double> y_coords_cache_;          ///< Cached Y coordinates
    mutable bool x_coords_cached_;                         ///< X coordinates cached flag
    mutable bool y_coords_cached_;                         ///< Y coordinates cached flag

    /**
     * @brief Invalidate coordinate cache
     */
    void invalidate_cache() noexcept {
        x_coords_cached_ = false;
        y_coords_cached_ = false;
    }

    /**
     * @brief Compute X coordinates and cache them
     */
    void compute_x_coords() const;

    /**
     * @brief Compute Y coordinates and cache them
     */
    void compute_y_coords() const;

    /**
     * @brief Validate grid parameters
     * @throws std::invalid_argument if parameters are invalid
     */
    void validate_parameters() const;
};

} // namespace mesh

#endif // COORDINATE_SYSTEM_HPP
