#include "domain_decomposition.hpp"
#include <iostream>
#include <sstream>
#include <algorithm>

// ========== Constructors ==========

DomainDecomposition::DomainDecomposition(size_t global_nx, size_t global_ny,
                                       const CartesianTopology& topology)
    : global_nx_(global_nx),
      global_ny_(global_ny),
      dim_x_(0),
      dim_y_(0),
      my_rank_(static_cast<size_t>(topology.rank())),
      topology_(&topology) {

    validate_topology();

    // Compute optimal dimensions
    auto dims = compute_optimal_dimensions(static_cast<size_t>(topology.size()),
                                           global_nx, global_ny);
    dim_x_ = dims.first;
    dim_y_ = dims.second;

    // Create Cartesian decomposition
    create_cartesian_decomposition();
}

DomainDecomposition::DomainDecomposition(size_t global_nx, size_t global_ny,
                                       size_t dim_x, size_t dim_y,
                                       const CartesianTopology& topology)
    : global_nx_(global_nx),
      global_ny_(global_ny),
      dim_x_(dim_x),
      dim_y_(dim_y),
      my_rank_(static_cast<size_t>(topology.rank())),
      topology_(&topology) {

    validate_topology();

    // Validate dimensions
    size_t dim_product = dim_x * dim_y;
    size_t num_procs = static_cast<size_t>(topology.size());
    if (dim_product != num_procs) {
        std::ostringstream oss;
        oss << "Invalid dimensions: product " << dim_x << " * " << dim_y << " = "
            << dim_product << " does not equal number of processes " << num_procs;
        throw std::invalid_argument(oss.str());
    }

    // Create Cartesian decomposition
    create_cartesian_decomposition();
}

DomainDecomposition::DomainDecomposition(std::vector<Subdomain> subdomains,
                                       const CartesianTopology& topology)
    : global_nx_(0),
      global_ny_(0),
      dim_x_(static_cast<size_t>(topology.dims()[0])),
      dim_y_(static_cast<size_t>(topology.dims()[1])),
      my_rank_(static_cast<size_t>(topology.rank())),
      subdomains_(std::move(subdomains)),
      topology_(&topology) {

    validate_topology();

    // Validate number of subdomains
    size_t num_procs = static_cast<size_t>(topology.size());
    if (subdomains_.size() != num_procs) {
        std::ostringstream oss;
        oss << "Invalid number of subdomains: expected " << num_procs
            << ", got " << subdomains_.size();
        throw std::invalid_argument(oss.str());
    }

    // Validate manual decomposition
    validate_manual_decomposition();

    // Set global dimensions from first subdomain
    if (!subdomains_.empty()) {
        global_nx_ = subdomains_[0].global_nx;
        global_ny_ = subdomains_[0].global_ny;
    }
}

// ========== Decomposition Query ==========

bool DomainDecomposition::is_in_my_domain(size_t global_i, size_t global_j) const {
    const auto& my_domain = my_subdomain();
    return (global_i >= my_domain.start_i &&
            global_i < my_domain.start_i + my_domain.nx &&
            global_j >= my_domain.start_j &&
            global_j < my_domain.start_j + my_domain.ny);
}

int DomainDecomposition::find_owner_rank(size_t global_i, size_t global_j) const {
    // Check bounds
    if (global_i >= global_nx_ || global_j >= global_ny_) {
        std::ostringstream oss;
        oss << "Global coordinates (" << global_i << ", " << global_j
            << ") are out of bounds (global_nx=" << global_nx_
            << ", global_ny=" << global_ny_ << ")";
        throw std::out_of_range(oss.str());
    }

    // For Cartesian decomposition, we can compute the owner directly
    // Find which subdomain contains these coordinates
    for (size_t rank = 0; rank < subdomains_.size(); ++rank) {
        const auto& subdomain = subdomains_[rank];
        if (global_i >= subdomain.start_i &&
            global_i < subdomain.start_i + subdomain.nx &&
            global_j >= subdomain.start_j &&
            global_j < subdomain.start_j + subdomain.ny) {
            return static_cast<int>(rank);
        }
    }

    // Should never reach here if decomposition is valid
    throw std::runtime_error("Failed to find owner for coordinates - invalid decomposition");
}

bool DomainDecomposition::is_on_boundary(Direction dir) const {
    validate_topology();

    switch (dir) {
        case Direction::X:
            // Check if on west or east boundary
            return topology_->is_on_boundary(Direction::X, Shift::Backward) ||
                   topology_->is_on_boundary(Direction::X, Shift::Forward);
        case Direction::Y:
            // Check if on south or north boundary
            return topology_->is_on_boundary(Direction::Y, Shift::Backward) ||
                   topology_->is_on_boundary(Direction::Y, Shift::Forward);
    }
    return false;
}

bool DomainDecomposition::has_neighbor(Direction dir) const {
    return !is_on_boundary(dir);
}

// ========== Load Analysis ==========

LoadBalanceStats DomainDecomposition::compute_load_balance() const {
    LoadBalanceStats stats;

    if (subdomains_.empty()) {
        return stats;
    }

    // Compute load for each subdomain
    std::vector<double> loads;
    loads.reserve(subdomains_.size());

    double total_load = 0.0;
    for (const auto& subdomain : subdomains_) {
        double load = static_cast<double>(subdomain.nx) * static_cast<double>(subdomain.ny);
        loads.push_back(load);
        total_load += load;
    }

    // Compute statistics
    stats.min_load = *std::min_element(loads.begin(), loads.end());
    stats.max_load = *std::max_element(loads.begin(), loads.end());
    stats.avg_load = total_load / static_cast<double>(loads.size());

    // Compute imbalance ratio
    if (stats.avg_load > 0.0) {
        stats.imbalance_ratio = (stats.max_load - stats.min_load) / stats.avg_load;
    }

    return stats;
}

// ========== Utility Methods ==========

std::pair<size_t, size_t> DomainDecomposition::compute_optimal_dimensions(
    size_t num_procs, size_t global_nx, size_t global_ny) {

    if (num_procs == 0) {
        return {0, 0};
    }

    if (num_procs == 1) {
        return {1, 1};
    }

    // Find dimensions that minimize aspect ratio difference
    // and balance load as evenly as possible
    size_t best_dim_x = 1;
    size_t best_dim_y = num_procs;
    double best_score = std::numeric_limits<double>::max();

    for (size_t dim_x = 1; dim_x <= num_procs; ++dim_x) {
        if (num_procs % dim_x == 0) {
            size_t dim_y = num_procs / dim_x;

            // Calculate the local domain sizes
            size_t local_nx = (global_nx + dim_x - 1) / dim_x;
            size_t local_ny = (global_ny + dim_y - 1) / dim_y;

            // Score based on aspect ratio (prefer square local domains)
            double aspect_ratio = static_cast<double>(local_nx) / static_cast<double>(local_ny);
            double score = std::abs(aspect_ratio - 1.0);

            if (score < best_score) {
                best_score = score;
                best_dim_x = dim_x;
                best_dim_y = dim_y;
            }
        }
    }

    return {best_dim_x, best_dim_y};
}

// ========== Private Methods ==========

void DomainDecomposition::create_cartesian_decomposition() {
    validate_topology();

    subdomains_.clear();
    subdomains_.reserve(dim_x_ * dim_y_);

    // Calculate base size of each subdomain
    size_t base_nx = global_nx_ / dim_x_;
    size_t base_ny = global_ny_ / dim_y_;

    // Calculate remainder (extra cells to distribute)
    size_t rem_nx = global_nx_ % dim_x_;
    size_t rem_ny = global_ny_ % dim_y_;

    // Create subdomains for each process
    size_t start_j = 0;
    for (size_t j = 0; j < dim_y_; ++j) {
        // Calculate Y size for this row (first rem_ny rows get +1)
        size_t ny = base_ny + (j < rem_ny ? 1 : 0);
        size_t start_i = 0;

        for (size_t i = 0; i < dim_x_; ++i) {
            // Calculate X size for this column (first rem_nx columns get +1)
            size_t nx = base_nx + (i < rem_nx ? 1 : 0);

            Subdomain subdomain;
            subdomain.start_i = start_i;
            subdomain.start_j = start_j;
            subdomain.nx = nx;
            subdomain.ny = ny;
            subdomain.global_nx = global_nx_;
            subdomain.global_ny = global_ny_;

            subdomains_.push_back(subdomain);

            start_i += nx;
        }

        start_j += ny;
    }
}

void DomainDecomposition::validate_manual_decomposition() const {
    if (subdomains_.empty()) {
        return;
    }

    // Create a 2D array to track which cells are assigned
    std::vector<std::vector<bool>> assigned(global_ny_, std::vector<bool>(global_nx_, false));

    // Check each subdomain
    for (const auto& subdomain : subdomains_) {
        // Validate global dimensions match
        if (subdomain.global_nx != global_nx_ || subdomain.global_ny != global_ny_) {
            std::ostringstream oss;
            oss << "Subdomain has invalid global dimensions: ("
                << subdomain.global_nx << ", " << subdomain.global_ny
                << ") != (" << global_nx_ << ", " << global_ny_ << ")";
            throw std::invalid_argument(oss.str());
        }

        // Validate subdomain is within bounds
        if (subdomain.start_i + subdomain.nx > global_nx_ ||
            subdomain.start_j + subdomain.ny > global_ny_) {
            std::ostringstream oss;
            oss << "Subdomain extends beyond global domain: start_i="
                << subdomain.start_i << ", nx=" << subdomain.nx
                << ", start_j=" << subdomain.start_j << ", ny=" << subdomain.ny;
            throw std::invalid_argument(oss.str());
        }

        // Mark cells as assigned and check for overlap
        for (size_t j = subdomain.start_j; j < subdomain.start_j + subdomain.ny; ++j) {
            for (size_t i = subdomain.start_i; i < subdomain.start_i + subdomain.nx; ++i) {
                if (assigned[j][i]) {
                    std::ostringstream oss;
                    oss << "Subdomain overlap detected at cell (" << i << ", " << j << ")";
                    throw std::invalid_argument(oss.str());
                }
                assigned[j][i] = true;
            }
        }
    }

    // Check that all cells are assigned
    for (size_t j = 0; j < global_ny_; ++j) {
        for (size_t i = 0; i < global_nx_; ++i) {
            if (!assigned[j][i]) {
                std::ostringstream oss;
                oss << "Cell (" << i << ", " << j << ") is not assigned to any subdomain";
                throw std::invalid_argument(oss.str());
            }
        }
    }
}

void DomainDecomposition::validate_topology() const {
    if (topology_ == nullptr) {
        throw std::runtime_error("Cartesian topology is null");
    }

    if (topology_->communicator() == MPI_COMM_NULL) {
        throw std::runtime_error("Cartesian topology has invalid communicator");
    }
}
