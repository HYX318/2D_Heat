/**
 * @file mesh2d.cpp
 * @brief Implementation of Mesh2D class
 */

#include "mesh2d.hpp"
#include <algorithm>
#include <stdexcept>

// Basic constructor - creates local mesh without ghost cells
Mesh2D::Mesh2D(size_t nx, size_t ny, double lx, double ly)
    : nx_(nx), ny_(ny),
      total_nx_(nx), total_ny_(ny),
      lx_(lx), ly_(ly),
      hx_(0.0), hy_(0.0),
      has_ghost_cells_(false),
      data_(nx, ny),
      topology_(nullptr),
      ghost_exchange_(nullptr),
      offset_x_(0), offset_y_(0) {

    if (nx == 0 || ny == 0) {
        throw std::invalid_argument("Mesh dimensions must be non-zero");
    }

    compute_spacings();
}

// Constructor with MPI - creates mesh with ghost cells
Mesh2D::Mesh2D(size_t nx, size_t ny, double lx, double ly,
               const CartesianTopology& topology)
    : nx_(nx), ny_(ny),
      total_nx_(nx + 2), total_ny_(ny + 2),
      lx_(lx), ly_(ly),
      hx_(0.0), hy_(0.0),
      has_ghost_cells_(true),
      data_(ny + 2, nx + 2),
      topology_(&topology),
      ghost_exchange_(nullptr),
      offset_x_(1), offset_y_(1) {

    if (nx == 0 || ny == 0) {
        throw std::invalid_argument("Mesh dimensions must be non-zero");
    }

    compute_spacings();
    initialize_ghost_exchange(topology);
}

// Copy constructor
Mesh2D::Mesh2D(const Mesh2D& other)
    : nx_(other.nx_), ny_(other.ny_),
      total_nx_(other.total_nx_), total_ny_(other.total_ny_),
      lx_(other.lx_), ly_(other.ly_),
      hx_(other.hx_), hy_(other.hy_),
      has_ghost_cells_(other.has_ghost_cells_),
      data_(other.data_),
      topology_(other.topology_),
      ghost_exchange_(nullptr),  // Don't copy ghost exchange
      offset_x_(other.offset_x_), offset_y_(other.offset_y_) {

    // If original had ghost exchange, create a new one
    if (other.ghost_exchange_ && topology_) {
        initialize_ghost_exchange(*topology_);
    }
}

// Move constructor
Mesh2D::Mesh2D(Mesh2D&& other) noexcept
    : nx_(other.nx_), ny_(other.ny_),
      total_nx_(other.total_nx_), total_ny_(other.total_ny_),
      lx_(other.lx_), ly_(other.ly_),
      hx_(other.hx_), hy_(other.hy_),
      has_ghost_cells_(other.has_ghost_cells_),
      data_(std::move(other.data_)),
      topology_(other.topology_),
      ghost_exchange_(std::move(other.ghost_exchange_)),
      offset_x_(other.offset_x_), offset_y_(other.offset_y_) {

    other.nx_ = 0;
    other.ny_ = 0;
    other.total_nx_ = 0;
    other.total_ny_ = 0;
    other.lx_ = 0.0;
    other.ly_ = 0.0;
    other.hx_ = 0.0;
    other.hy_ = 0.0;
    other.has_ghost_cells_ = false;
    other.topology_ = nullptr;
    other.offset_x_ = 0;
    other.offset_y_ = 0;
}

// Destructor
Mesh2D::~Mesh2D() = default;

// Copy assignment operator
Mesh2D& Mesh2D::operator=(const Mesh2D& other) {
    if (this != &other) {
        Mesh2D temp(other);
        std::swap(nx_, temp.nx_);
        std::swap(ny_, temp.ny_);
        std::swap(total_nx_, temp.total_nx_);
        std::swap(total_ny_, temp.total_ny_);
        std::swap(lx_, temp.lx_);
        std::swap(ly_, temp.ly_);
        std::swap(hx_, temp.hx_);
        std::swap(hy_, temp.hy_);
        std::swap(has_ghost_cells_, temp.has_ghost_cells_);
        std::swap(data_, temp.data_);
        std::swap(topology_, temp.topology_);
        std::swap(ghost_exchange_, temp.ghost_exchange_);
        std::swap(offset_x_, temp.offset_x_);
        std::swap(offset_y_, temp.offset_y_);
    }
    return *this;
}

// Move assignment operator
Mesh2D& Mesh2D::operator=(Mesh2D&& other) noexcept {
    if (this != &other) {
        nx_ = other.nx_;
        ny_ = other.ny_;
        total_nx_ = other.total_nx_;
        total_ny_ = other.total_ny_;
        lx_ = other.lx_;
        ly_ = other.ly_;
        hx_ = other.hx_;
        hy_ = other.hy_;
        has_ghost_cells_ = other.has_ghost_cells_;
        data_ = std::move(other.data_);
        topology_ = other.topology_;
        ghost_exchange_ = std::move(other.ghost_exchange_);
        offset_x_ = other.offset_x_;
        offset_y_ = other.offset_y_;

        other.nx_ = 0;
        other.ny_ = 0;
        other.total_nx_ = 0;
        other.total_ny_ = 0;
        other.lx_ = 0.0;
        other.ly_ = 0.0;
        other.hx_ = 0.0;
        other.hy_ = 0.0;
        other.has_ghost_cells_ = false;
        other.topology_ = nullptr;
        other.offset_x_ = 0;
        other.offset_y_ = 0;
    }
    return *this;
}

// Element access (mutable) - accounts for ghost cell offset
double& Mesh2D::operator()(size_t i, size_t j) {
    check_interior_bounds(i, j);
    return data_(offset_y_ + i, offset_x_ + j);
}

// Element access (const) - accounts for ghost cell offset
const double& Mesh2D::operator()(size_t i, size_t j) const {
    check_interior_bounds(i, j);
    return data_(offset_y_ + i, offset_x_ + j);
}

// Direct element access to underlying array (mutable)
double& Mesh2D::at(size_t i, size_t j) {
    check_total_bounds(i, j);
    return data_(i, j);
}

// Direct element access to underlying array (const)
const double& Mesh2D::at(size_t i, size_t j) const {
    check_total_bounds(i, j);
    return data_(i, j);
}

// Convert global x-index to local x-index
size_t Mesh2D::global_to_local_x(size_t global_i) const {
    if (!topology_) {
        // No MPI - global = local
        return global_i;
    }

    // Calculate local index based on topology
    int coord_x = topology_->coord_x();
    int dim_x = topology_->dim_x();
    size_t local_nx = nx_ / dim_x;

    if (global_i < coord_x * local_nx || global_i >= (coord_x + 1) * local_nx) {
        throw std::out_of_range("Global index is not in local domain");
    }

    return global_i - coord_x * local_nx;
}

// Convert global y-index to local y-index
size_t Mesh2D::global_to_local_y(size_t global_j) const {
    if (!topology_) {
        // No MPI - global = local
        return global_j;
    }

    // Calculate local index based on topology
    int coord_y = topology_->coord_y();
    int dim_y = topology_->dim_y();
    size_t local_ny = ny_ / dim_y;

    if (global_j < coord_y * local_ny || global_j >= (coord_y + 1) * local_ny) {
        throw std::out_of_range("Global index is not in local domain");
    }

    return global_j - coord_y * local_ny;
}

// Convert local x-index to global x-index
size_t Mesh2D::local_to_global_x(size_t local_i) const {
    if (!topology_) {
        return local_i;
    }

    int coord_x = topology_->coord_x();
    int dim_x = topology_->dim_x();
    size_t local_nx = nx_ / dim_x;

    return coord_x * local_nx + local_i;
}

// Convert local y-index to global y-index
size_t Mesh2D::local_to_global_y(size_t local_j) const {
    if (!topology_) {
        return local_j;
    }

    int coord_y = topology_->coord_y();
    int dim_y = topology_->dim_y();
    size_t local_ny = ny_ / dim_y;

    return coord_y * local_ny + local_j;
}

// Get physical x-coordinate at grid point
double Mesh2D::x_coord(size_t i) const {
    check_x_bounds(i);
    return static_cast<double>(i) * hx_;
}

// Get physical y-coordinate at grid point
double Mesh2D::y_coord(size_t j) const {
    check_y_bounds(j);
    return static_cast<double>(j) * hy_;
}

// Get physical coordinates at grid point
std::pair<double, double> Mesh2D::coord(size_t i, size_t j) const {
    return std::make_pair(x_coord(i), y_coord(j));
}

// Apply constant Dirichlet boundary condition
void Mesh2D::apply_dirichlet_bc(double value, double t) {
    // Time parameter is ignored for constant BC
    (void)t;  // Suppress unused parameter warning

    if (!has_ghost_cells_) {
        // No ghost cells - no interior points to apply BC to
        // For non-ghost mesh, we might want to apply BC to boundary indices
        // but that depends on the problem setup
        return;
    }

    // Set all ghost cells to the boundary value
    set_boundary_values(value);
}

// Apply function-based boundary condition
void Mesh2D::apply_bc(std::function<double(double x, double y, double t)> bc_func, double t) {
    if (!has_ghost_cells_) {
        return;
    }

    // Apply to all ghost cell boundaries
    // South ghost row
    for (size_t j = 0; j < total_nx_; ++j) {
        double x = j * hx_;
        double y = 0.0;
        data_(0, j) = bc_func(x, y, t);
    }

    // North (ny+1) ghost row
    for (size_t j = 0; j < total_nx_; ++j) {
        double x = j * hx_;
        double y = ly_;
        data_(total_ny_ - 1, j) = bc_func(x, y, t);
    }

    // West ghost column
    for (size_t i = 0; i < total_ny_; ++i) {
        double x = 0.0;
        double y = i * hy_;
        data_(i, 0) = bc_func(x, y, t);
    }

    // East (nx+1) ghost column
    for (size_t i = 0; i < total_ny_; ++i) {
        double x = lx_;
        double y = i * hy_;
        data_(i, total_nx_ - 1) = bc_func(x, y, t);
    }
}

// Apply constant boundary condition at specific direction
void Mesh2D::apply_bc_at_direction(CardinalDirection dir, double value) {
    if (!has_ghost_cells_) {
        return;
    }

    switch (dir) {
        case CardinalDirection::South:
            // South ghost row
            for (size_t j = 0; j < total_nx_; ++j) {
                data_(0, j) = value;
            }
            break;

        case CardinalDirection::North:
            // North (ny+1) ghost row
            for (size_t j = 0; j < total_nx_; ++j) {
                data_(total_ny_ - 1, j) = value;
            }
            break;

        case CardinalDirection::West:
            // West ghost column
            for (size_t i = 0; i < total_ny_; ++i) {
                data_(i, 0) = value;
            }
            break;

        case CardinalDirection::East:
            // East (nx+1) ghost column
            for (size_t i = 0; i < total_ny_; ++i) {
                data_(i, total_nx_ - 1) = value;
            }
            break;
    }
}

// Compute Laplacian using 5-point stencil
void Mesh2D::compute_laplacian(utils::Array2D& result) const {
    if (result.rows() != ny_ || result.cols() != nx_) {
        throw std::invalid_argument("Result array dimensions must match mesh interior dimensions");
    }

    double inv_hx2 = 1.0 / (hx_ * hx_);
    double inv_hy2 = 1.0 / (hy_ * hy_);

    for (size_t i = 0; i < ny_; ++i) {
        for (size_t j = 0; j < nx_; ++j) {
            // Convert to array indices (accounting for ghost cells)
            size_t ai = offset_y_ + i;
            size_t aj = offset_x_ + j;

            // Central difference in x
            double d2x = (data_(ai, aj + 1) - 2.0 * data_(ai, aj) + data_(ai, aj - 1)) * inv_hx2;

            // Central difference in y
            double d2y = (data_(ai + 1, aj) - 2.0 * data_(ai, aj) + data_(ai - 1, aj)) * inv_hy2;

            result(i, j) = d2x + d2y;
        }
    }
}

// Fill all interior points with specified value
void Mesh2D::fill(double value) {
    for (size_t i = 0; i < ny_; ++i) {
        for (size_t j = 0; j < nx_; ++j) {
            (*this)(i, j) = value;
        }
    }
}

// Copy data from another mesh
void Mesh2D::copy_from(const Mesh2D& other) {
    if (nx_ != other.nx_ || ny_ != other.ny_) {
        throw std::invalid_argument("Mesh dimensions must match for copy_from");
    }

    for (size_t i = 0; i < ny_; ++i) {
        for (size_t j = 0; j < nx_; ++j) {
            (*this)(i, j) = other(i, j);
        }
    }
}

// Scale all interior points by factor
void Mesh2D::scale(double factor) {
    for (size_t i = 0; i < ny_; ++i) {
        for (size_t j = 0; j < nx_; ++j) {
            (*this)(i, j) *= factor;
        }
    }
}

// Add another mesh to this mesh
void Mesh2D::add(const Mesh2D& other) {
    if (nx_ != other.nx_ || ny_ != other.ny_) {
        throw std::invalid_argument("Mesh dimensions must match for addition");
    }

    for (size_t i = 0; i < ny_; ++i) {
        for (size_t j = 0; j < nx_; ++j) {
            (*this)(i, j) += other(i, j);
        }
    }
}

// Subtract another mesh from this mesh
void Mesh2D::subtract(const Mesh2D& other) {
    if (nx_ != other.nx_ || ny_ != other.ny_) {
        throw std::invalid_argument("Mesh dimensions must match for subtraction");
    }

    for (size_t i = 0; i < ny_; ++i) {
        for (size_t j = 0; j < nx_; ++j) {
            (*this)(i, j) -= other(i, j);
        }
    }
}

// Compute L2 norm (Euclidean norm)
double Mesh2D::l2_norm() const {
    double sum = 0.0;
    for (size_t i = 0; i < ny_; ++i) {
        for (size_t j = 0; j < nx_;
             ++j) {
            double val = (*this)(i, j);
            sum += val * val;
        }
    }
    return std::sqrt(sum);
}

// Compute L-infinity norm (maximum absolute value)
double Mesh2D::linfty_norm() const {
    double max_abs = 0.0;
    for (size_t i = 0; i < ny_; ++i) {
        for (size_t j = 0; j < nx_; ++j) {
            max_abs = std::max(max_abs, std::abs((*this)(i, j)));
        }
    }
    return max_abs;
}

// Get maximum value
double Mesh2D::max() const {
    double max_val = (*this)(0, 0);
    for (size_t i = 0; i < ny_; ++i) {
        for (size_t j = 0; j < nx_; ++j) {
            max_val = std::max(max_val, (*this)(i, j));
        }
    }
    return max_val;
}

// Get minimum value
double Mesh2D::min() const {
    double min_val = (*this)(0, 0);
    for (size_t i = 0; i < ny_; ++i) {
        for (size_t j = 0; j < nx_; ++j) {
            min_val = std::min(min_val, (*this)(i, j));
        }
    }
    return min_val;
}

// Exchange ghost cells with neighboring processes
void Mesh2D::exchange_ghost_cells() {
    if (!has_ghost_cells_) {
        throw std::runtime_error("Ghost cells are not enabled for this mesh");
    }

    if (!ghost_exchange_) {
        throw std::runtime_error("Ghost cell exchange is not initialized");
    }

    ghost_exchange_->exchange(data_);
}

// Set all boundary values (including ghost cells)
void Mesh2D::set_boundary_values(double value) {
    if (!has_ghost_cells_) {
        return;
    }

    // South ghost row
    for (size_t j = 0; j < total_nx_; ++j) {
        data_(0, j) = value;
    }

    // North (ny+1) ghost row
    for (size_t j = 0; j < total_nx_; ++j) {
        data_(total_ny_ - 1, j) = value;
    }

    // West ghost column
    for (size_t i = 0; i < total_ny_; ++i) {
        data_(i, 0) = value;
    }

    // East (nx+1) ghost column
    for (size_t i = 0; i < total_ny_; ++i) {
        data_(i, total_nx_ - 1) = value;
    }
}

// Validate indices are within interior bounds
void Mesh2D::check_interior_bounds(size_t i, size_t j) const {
    if (i >= ny_) {
        throw std::out_of_range("Row index out of interior bounds");
    }
    if (j >= nx_) {
        throw std::out_of_range("Column index out of interior bounds");
    }
}

// Validate indices are within total bounds (including ghost)
void Mesh2D::check_total_bounds(size_t i, size_t j) const {
    if (i >= total_ny_) {
        throw std::out_of_range("Row index out of total bounds");
    }
    if (j >= total_nx_) {
        throw std::out_of_range("Column index out of total bounds");
    }
}

// Validate index is within x bounds
void Mesh2D::check_x_bounds(size_t i) const {
    if (i >= nx_) {
        throw std::out_of_range("X-index out of bounds");
    }
}

// Validate index is within y bounds
void Mesh2D::check_y_bounds(size_t j) const {
    if (j >= ny_) {
        throw std::out_of_range("Y-index out of bounds");
    }
}

// Compute grid spacings
void Mesh2D::compute_spacings() {
    // Grid spacing: hx = lx / (nx - 1) for inclusive endpoints
    // or hx = lx / nx for cell-centered
    // We'll use inclusive endpoint interpretation
    if (nx_ > 1) {
        hx_ = lx_ / static_cast<double>(nx_ - 1);
    } else {
        hx_ = lx_;
    }

    if (ny_ > 1) {
        hy_ = ly_ / static_cast<double>(ny_ - 1);
    } else {
        hy_ = ly_;
    }
}

// Initialize ghost cell exchange
void Mesh2D::initialize_ghost_exchange(const CartesianTopology& topology) {
    try {
        ghost_exchange_ = std::make_unique<GhostCellExchange>(
            static_cast<int>(nx_),
            static_cast<int>(ny_),
            topology
        );
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("Failed to initialize ghost cell exchange: ") + e.what());
    }
}
