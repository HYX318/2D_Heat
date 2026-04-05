#include "cartesian_topology.hpp"
#include <iostream>
#include <sstream>
#include <algorithm>

CartesianTopology::CartesianTopology(MPI_Comm comm)
    : cart_comm_(MPI_COMM_NULL),
      parent_comm_(comm),
      rank_(-1),
      size_(0),
      owns_comm_(true) {

    // Get the number of processes in the parent communicator
    int mpi_comm_size_result = MPI_Comm_size(parent_comm_, &size_);
    if (mpi_comm_size_result != MPI_SUCCESS) {
        std::ostringstream oss;
        oss << "MPI_Comm_size failed with error code: " << mpi_comm_size_result;
        throw std::runtime_error(oss.str());
    }

    // Compute optimal dimensions
    dims_ = compute_optimal_dimensions(size_);

    // Create periods array (non-periodic by default)
    std::vector<int> periods(2, 0);

    // Reorder is true (allow MPI to reassign ranks for performance)
    int reorder = 1;

    // Create the Cartesian communicator
    int mpi_cart_create_result = MPI_Cart_create(
        parent_comm_, 2, dims_.data(), periods.data(), reorder, &cart_comm_
    );

    if (mpi_cart_create_result != MPI_SUCCESS) {
        std::ostringstream oss;
        oss << "MPI_Cart_create failed with error code: " << mpi_cart_create_result;
        throw std::runtime_error(oss.str());
    }

    // Get process rank
    int mpi_comm_rank_result = MPI_Comm_rank(cart_comm_, &rank_);
    if (mpi_comm_rank_result != MPI_SUCCESS) {
        std::ostringstream oss;
        oss << "MPI_Comm_rank failed with error code: " << mpi_comm_rank_result;
        throw std::runtime_error(oss.str());
    }

    // Get process coordinates
    coords_.resize(2);
    int mpi_cart_coords_result = MPI_Cart_coords(cart_comm_, rank_, 2, coords_.data());
    if (mpi_cart_coords_result != MPI_SUCCESS) {
        std::ostringstream oss;
        oss << "MPI_Cart_coords failed with error code: " << mpi_cart_coords_result;
        throw std::runtime_error(oss.str());
    }

    // Compute neighbors
    compute_neighbors();
}

CartesianTopology::CartesianTopology(MPI_Comm comm, const std::vector<int>& dims)
    : cart_comm_(MPI_COMM_NULL),
      parent_comm_(comm),
      rank_(-1),
      size_(0),
      owns_comm_(true) {

    // Validate dimensions
    if (dims.size() != 2) {
        std::ostringstream oss;
        oss << "Invalid dimensions: expected 2 elements, got " << dims.size();
        throw std::invalid_argument(oss.str());
    }

    // Get the number of processes in the parent communicator
    int mpi_comm_size_result = MPI_Comm_size(parent_comm_, &size_);
    if (mpi_comm_size_result != MPI_SUCCESS) {
        std::ostringstream oss;
        oss << "MPI_Comm_size failed with error code: " << mpi_comm_size_result;
        throw std::runtime_error(oss.str());
    }

    // Validate that product of dimensions equals number of processes
    int dim_product = dims[0] * dims[1];
    if (dim_product != size_) {
        std::ostringstream oss;
        oss << "Invalid dimensions: product " << dims[0] << " * " << dims[1] << " = "
            << dim_product << " does not equal number of processes " << size_;
        throw std::invalid_argument(oss.str());
    }

    dims_ = dims;

    // Create periods array (non-periodic by default)
    std::vector<int> periods(2, 0);

    // Reorder is true (allow MPI to reassign ranks for performance)
    int reorder = 1;

    // Create the Cartesian communicator
    int mpi_cart_create_result = MPI_Cart_create(
        parent_comm_, 2, dims_.data(), periods.data(), reorder, &cart_comm_
    );

    if (mpi_cart_create_result != MPI_SUCCESS) {
        std::ostringstream oss;
        oss << "MPI_Cart_create failed with error code: " << mpi_cart_create_result;
        throw std::runtime_error(oss.str());
    }

    // Get process rank
    int mpi_comm_rank_result = MPI_Comm_rank(cart_comm_, &rank_);
    if (mpi_comm_rank_result != MPI_SUCCESS) {
        std::ostringstream oss;
        oss << "MPI_Comm_rank failed with error code: " << mpi_comm_rank_result;
        throw std::runtime_error(oss.str());
    }

    // Get process coordinates
    coords_.resize(2);
    int mpi_cart_coords_result = MPI_Cart_coords(cart_comm_, rank_, 2, coords_.data());
    if (mpi_cart_coords_result != MPI_SUCCESS) {
        std::ostringstream oss;
        oss << "MPI_Cart_coords failed with error code: " << mpi_cart_coords_result;
        throw std::runtime_error(oss.str());
    }

    // Compute neighbors
    compute_neighbors();
}

CartesianTopology::~CartesianTopology() {
    try {
        // Free the Cartesian communicator if we own it
        if (cart_comm_ != MPI_COMM_NULL && cart_comm_ != MPI_COMM_WORLD && owns_comm_) {
            int mpi_comm_free_result = MPI_Comm_free(&cart_comm_);
            if (mpi_comm_free_result != MPI_SUCCESS) {
                std::cerr << "Warning: MPI_Comm_free failed with error code: "
                          << mpi_comm_free_result << std::endl;
            }
            cart_comm_ = MPI_COMM_NULL;
        }
    } catch (...) {
        // Destructor should never throw - catch any exceptions
        std::cerr << "Warning: Exception caught in CartesianTopology destructor" << std::endl;
    }
}

CartesianTopology::CartesianTopology(CartesianTopology&& other) noexcept
    : cart_comm_(other.cart_comm_),
      parent_comm_(other.parent_comm_),
      rank_(other.rank_),
      size_(other.size_),
      dims_(std::move(other.dims_)),
      coords_(std::move(other.coords_)),
      neighbors_(other.neighbors_),
      owns_comm_(other.owns_comm_) {

    // Mark the moved-from object as invalid
    other.cart_comm_ = MPI_COMM_NULL;
    other.owns_comm_ = false;
    other.rank_ = -1;
    other.size_ = 0;
    other.dims_.clear();
    other.coords_.clear();
    other.neighbors_ = NeighborInfo();
}

CartesianTopology& CartesianTopology::operator=(CartesianTopology&& other) noexcept {
    if (this != &other) {
        // Clean up existing resources if we own them
        if (cart_comm_ != MPI_COMM_NULL && cart_comm_ != MPI_COMM_WORLD && owns_comm_) {
            try {
                int mpi_comm_free_result = MPI_Comm_free(&cart_comm_);
                if (mpi_comm_free_result != MPI_SUCCESS) {
                    std::cerr << "Warning: MPI_Comm_free failed with error code: "
                              << mpi_comm_free_result << std::endl;
                }
            } catch (...) {
                // Assignment operator should never throw
                std::cerr << "Warning: Exception caught in CartesianTopology move assignment"
                          << std::endl;
            }
        }

        // Transfer ownership
        cart_comm_ = other.cart_comm_;
        parent_comm_ = other.parent_comm_;
        rank_ = other.rank_;
        size_ = other.size_;
        dims_ = std::move(other.dims_);
        coords_ = std::move(other.coords_);
        neighbors_ = other.neighbors_;
        owns_comm_ = other.owns_comm_;

        // Mark the moved-from object as invalid
        other.cart_comm_ = MPI_COMM_NULL;
        other.owns_comm_ = false;
        other.rank_ = -1;
        other.size_ = 0;
        other.dims_.clear();
        other.coords_.clear();
        other.neighbors_ = NeighborInfo();
    }
    return *this;
}

bool CartesianTopology::is_on_boundary(Direction dir, Shift shift) const {
    validate_communicator();

    switch (dir) {
        case Direction::X:
            if (shift == Shift::Backward) {
                return coords_[0] == 0;  // West boundary
            } else {  // Shift::Forward
                return coords_[0] == dims_[0] - 1;  // East boundary
            }
        case Direction::Y:
            if (shift == Shift::Backward) {
                return coords_[1] == 0;  // South boundary
            } else {  // Shift::Forward
                return coords_[1] == dims_[1] - 1;  // North boundary
            }
    }
    return false;
}

bool CartesianTopology::has_neighbor(Direction dir, Shift shift) const {
    return !is_on_boundary(dir, shift);
}

std::vector<int> CartesianTopology::compute_optimal_dimensions(int num_procs) {
    // Find the divisor closest to sqrt(num_procs)
    int sqrt_procs = static_cast<int>(std::sqrt(num_procs));
    int best_dim_x = 1;
    int best_dim_y = num_procs;
    int min_diff = std::abs(best_dim_x - best_dim_y);

    for (int dim_x = sqrt_procs; dim_x >= 1; --dim_x) {
        if (num_procs % dim_x == 0) {
            int dim_y = num_procs / dim_x;
            int diff = std::abs(dim_x - dim_y);
            if (diff < min_diff) {
                min_diff = diff;
                best_dim_x = dim_x;
                best_dim_y = dim_y;
            }
            // Found the closest, can break
            break;
        }
    }

    return {best_dim_x, best_dim_y};
}

void CartesianTopology::compute_neighbors() {
    validate_communicator();

    int source, dest;

    // Compute north/south neighbors (Y direction)
    MPI_Cart_shift(cart_comm_, static_cast<int>(Direction::Y),
                   static_cast<int>(Shift::Forward), &source, &dest);
    neighbors_.north = dest;

    MPI_Cart_shift(cart_comm_, static_cast<int>(Direction::Y),
                   static_cast<int>(Shift::Backward), &source, &dest);
    neighbors_.south = dest;

    // Compute east/west neighbors (X direction)
    MPI_Cart_shift(cart_comm_, static_cast<int>(Direction::X),
                   static_cast<int>(Shift::Forward), &source, &dest);
    neighbors_.east = dest;

    MPI_Cart_shift(cart_comm_, static_cast<int>(Direction::X),
                   static_cast<int>(Shift::Backward), &source, &dest);
    neighbors_.west = dest;
}

void CartesianTopology::validate_communicator() const {
    if (cart_comm_ == MPI_COMM_NULL) {
        throw std::runtime_error("Cartesian communicator is invalid (MPI_COMM_NULL)");
    }
}
