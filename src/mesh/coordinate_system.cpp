/**
 * @file coordinate_system.cpp
 * @brief Implementation of CoordinateSystem class
 */

#include "coordinate_system.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace mesh {

// ========== Constructors ==========

CoordinateSystem::CoordinateSystem(size_t nx, size_t ny,
                                   double x_min, double x_max,
                                   double y_min, double y_max)
    : nx_(nx), ny_(ny),
      x_min_(x_min), x_max_(x_max),
      y_min_(y_min), y_max_(y_max),
      grid_type_(GridType::Uniform),
      x_distribution_(nullptr),
      y_distribution_(nullptr),
      x_coords_cached_(false),
      y_coords_cached_(false) {
    validate_parameters();
}

CoordinateSystem::CoordinateSystem(const GridParams& params)
    : nx_(params.nx), ny_(params.ny),
      x_min_(params.x_min), x_max_(params.x_max),
      y_min_(params.y_min), y_max_(params.y_max),
      grid_type_(params.type),
      x_distribution_(params.x_distribution),
      y_distribution_(params.y_distribution),
      x_coords_cached_(false),
      y_coords_cached_(false) {
    validate_parameters();
}

// ========== Coordinate Calculation ==========

double CoordinateSystem::x(size_t i) const {
    if (i >= nx_) {
        throw std::out_of_range("X index out of range");
    }

    if (grid_type_ == GridType::Uniform) {
        return x_min_ + i * hx();
    } else {
        // Use cached coordinates for non-uniform grids
        if (!x_coords_cached_) {
            compute_x_coords();
        }
        return x_coords_cache_[i];
    }
}

double CoordinateSystem::y(size_t j) const {
    if (j >= ny_) {
        throw std::out_of_range("Y index out of range");
    }

    if (grid_type_ == GridType::Uniform) {
        return y_min_ + j * hy();
    } else {
        // Use cached coordinates for non-uniform grids
        if (!y_coords_cached_) {
            compute_y_coords();
        }
        return y_coords_cache_[j];
    }
}

size_t CoordinateSystem::i_at(double x) const {
    if (grid_type_ == GridType::Uniform) {
        double xi = (x - x_min_) / lx();
        size_t i = static_cast<size_t>(std::round(xi * (nx_ - 1)));
        return std::min(i, nx_ - 1);
    } else {
        // For non-uniform grids, use binary search
        if (!x_coords_cached_) {
            compute_x_coords();
        }
        auto it = std::lower_bound(x_coords_cache_.begin(),
                                   x_coords_cache_.end(), x);
        size_t i = std::distance(x_coords_cache_.begin(), it);

        // Find closest point
        if (i > 0 && i < nx_) {
            double dist_prev = std::abs(x_coords_cache_[i-1] - x);
            double dist_curr = std::abs(x_coords_cache_[i] - x);
            if (dist_prev < dist_curr) {
                return i - 1;
            }
        }
        return std::min(i, nx_ - 1);
    }
}

size_t CoordinateSystem::j_at(double y) const {
    if (grid_type_ == GridType::Uniform) {
        double eta = (y - y_min_) / ly();
        size_t j = static_cast<size_t>(std::round(eta * (ny_ - 1)));
        return std::min(j, ny_ - 1);
    } else {
        // For non-uniform grids, use binary search
        if (!y_coords_cached_) {
            compute_y_coords();
        }
        auto it = std::lower_bound(y_coords_cache_.begin(),
                                   y_coords_cache_.end(), y);
        size_t j = std::distance(y_coords_cache_.begin(), it);

        // Find closest point
        if (j > 0 && j < ny_) {
            double dist_prev = std::abs(y_coords_cache_[j-1] - y);
            double dist_curr = std::abs(y_coords_cache_[j] - y);
            if (dist_prev < dist_curr) {
                return j - 1;
            }
        }
        return std::min(j, ny_ - 1);
    }
}

std::pair<double, double> CoordinateSystem::coord(size_t i, size_t j) const {
    return std::make_pair(x(i), y(j));
}

double CoordinateSystem::hx_at(size_t i) const {
    if (i >= nx_ - 1) {
        if (nx_ > 1) {
            return x(nx_ - 1) - x(nx_ - 2);
        } else {
            return 0.0;
        }
    }
    return x(i + 1) - x(i);
}

double CoordinateSystem::hy_at(size_t j) const {
    if (j >= ny_ - 1) {
        if (ny_ > 1) {
            return y(ny_ - 1) - y(ny_ - 2);
        } else {
            return 0.0;
        }
    }
    return y(j + 1) - y(j);
}

// ========== Coordinate Transformations ==========

double CoordinateSystem::xi(size_t i) const {
    if (nx_ <= 1) return 0.0;
    return static_cast<double>(i) / (nx_ - 1);
}

double CoordinateSystem::eta(size_t j) const {
    if (ny_ <= 1) return 0.0;
    return static_cast<double>(j) / (ny_ - 1);
}

double CoordinateSystem::xi_from_x(double x) const {
    return (x - x_min_) / lx();
}

double CoordinateSystem::eta_from_y(double y) const {
    return (y - y_min_) / ly();
}

double CoordinateSystem::x_from_xi(double xi_normalized) const {
    if (grid_type_ == GridType::Uniform || !x_distribution_) {
        return x_min_ + xi_normalized * lx();
    } else {
        // Use distribution function
        double xi_clamped = std::max(0.0, std::min(1.0, xi_normalized));
        double mapped_xi = x_distribution_(xi_clamped);
        return x_min_ + mapped_xi * lx();
    }
}

double CoordinateSystem::y_from_eta(double eta) const {
    if (grid_type_ == GridType::Uniform || !y_distribution_) {
        return y_min_ + eta * ly();
    } else {
        // Use distribution function
        double eta_clamped = std::max(0.0, std::min(1.0, eta));
        double mapped_eta = y_distribution_(eta_clamped);
        return y_min_ + mapped_eta * ly();
    }
}

// ========== Grid Generation ==========

std::vector<double> CoordinateSystem::generate_x_coords() const {
    std::vector<double> coords(nx_);
    for (size_t i = 0; i < nx_; ++i) {
        coords[i] = x(i);
    }
    return coords;
}

std::vector<double> CoordinateSystem::generate_y_coords() const {
    std::vector<double> coords(ny_);
    for (size_t j = 0; j < ny_; ++j) {
        coords[j] = y(j);
    }
    return coords;
}

void CoordinateSystem::generate_coordinate_mesh(utils::Array2D& x_mesh,
                                                 utils::Array2D& y_mesh) const {
    if (x_mesh.rows() != ny_ || x_mesh.cols() != nx_) {
        throw std::invalid_argument("x_mesh dimensions do not match grid");
    }
    if (y_mesh.rows() != ny_ || y_mesh.cols() != nx_) {
        throw std::invalid_argument("y_mesh dimensions do not match grid");
    }

    // Pre-compute coordinate arrays
    std::vector<double> x_coords = generate_x_coords();
    std::vector<double> y_coords = generate_y_coords();

    // Fill coordinate meshes
    for (size_t j = 0; j < ny_; ++j) {
        for (size_t i = 0; i < nx_; ++i) {
            x_mesh(j, i) = x_coords[i];
            y_mesh(j, i) = y_coords[j];
        }
    }
}

// ========== Non-Uniform Grid Support ==========

void CoordinateSystem::set_x_distribution(std::function<double(double)> dist) {
    x_distribution_ = dist;
    grid_type_ = GridType::NonUniform;
    invalidate_cache();
}

void CoordinateSystem::set_y_distribution(std::function<double(double)> dist) {
    y_distribution_ = dist;
    grid_type_ = GridType::NonUniform;
    invalidate_cache();
}

double CoordinateSystem::stretching_function(double xi, double beta) {
    // Clamp xi to [0, 1]
    xi = std::max(0.0, std::min(1.0, xi));

    // Avoid division by zero
    if (std::abs(beta - 1.0) < 1e-10) {
        return xi;
    }

    // Stretching function: ξ_stretched = (1 - β^ξ) / (1 - β)
    return (1.0 - std::pow(beta, xi)) / (1.0 - beta);
}

// ========== Validation and Utilities ==========

bool CoordinateSystem::is_inside(double x, double y) const {
    return is_inside(x, y, 0.0);
}

bool CoordinateSystem::is_inside(double x, double y, double tolerance) const {
    return (x >= x_min_ - tolerance && x <= x_max_ + tolerance &&
            y >= y_min_ - tolerance && y <= y_max_ + tolerance);
}

double CoordinateSystem::grid_quality() const {
    if (grid_type_ == GridType::Uniform) {
        return 1.0;
    }

    // Compute spacing variation for non-uniform grids
    std::vector<double> x_spacings, y_spacings;

    for (size_t i = 0; i < nx_ - 1; ++i) {
        x_spacings.push_back(hx_at(i));
    }
    for (size_t j = 0; j < ny_ - 1; ++j) {
        y_spacings.push_back(hy_at(j));
    }

    if (x_spacings.empty() || y_spacings.empty()) {
        return 1.0;
    }

    // Compute mean and standard deviation of spacings
    auto compute_quality = [](const std::vector<double>& spacings) -> double {
        double mean = 0.0;
        for (double s : spacings) {
            mean += s;
        }
        mean /= spacings.size();

        double variance = 0.0;
        for (double s : spacings) {
            variance += (s - mean) * (s - mean);
        }
        variance /= spacings.size();

        double stddev = std::sqrt(variance);
        if (std::abs(mean) < 1e-10) {
            return 1.0;
        }
        return std::max(0.0, 1.0 - stddev / mean);
    };

    double x_quality = compute_quality(x_spacings);
    double y_quality = compute_quality(y_spacings);

    return (x_quality + y_quality) / 2.0;
}

double CoordinateSystem::aspect_ratio() const {
    if (std::abs(ly()) < 1e-10) {
        return std::numeric_limits<double>::infinity();
    }
    return lx() / ly();
}

// ========== Private Methods ==========

void CoordinateSystem::compute_x_coords() const {
    x_coords_cache_.resize(nx_);
    for (size_t i = 0; i < nx_; ++i) {
        double xi_norm = xi(i);
        if (x_distribution_) {
            x_coords_cache_[i] = x_from_xi(xi_norm);
        } else {
            x_coords_cache_[i] = x_min_ + xi_norm * lx();
        }
    }
    x_coords_cached_ = true;
}

void CoordinateSystem::compute_y_coords() const {
    y_coords_cache_.resize(ny_);
    for (size_t j = 0; j < ny_; ++j) {
        double eta_norm = eta(j);
        if (y_distribution_) {
            y_coords_cache_[j] = y_from_eta(eta_norm);
        } else {
            y_coords_cache_[j] = y_min_ + eta_norm * ly();
        }
    }
    y_coords_cached_ = true;
}

void CoordinateSystem::validate_parameters() const {
    if (nx_ == 0) {
        throw std::invalid_argument("nx must be greater than 0");
    }
    if (ny_ == 0) {
        throw std::invalid_argument("ny must be greater than 0");
    }
    if (x_max_ <= x_min_) {
        throw std::invalid_argument("x_max must be greater than x_min");
    }
    if (y_max_ <= y_min_) {
        throw std::invalid_argument("y_max must be greater than y_min");
    }
}

} // namespace mesh
