#include "ghost_cell_exchange.hpp"
#include <iostream>
#include <string>

GhostCellExchange::GhostCellExchange(int nx, int ny, const CartesianTopology& topology)
    : nx_(nx), ny_(ny), topology(topology), row_type_(MPI_DATATYPE_NULL),
      col_type_(MPI_DATATYPE_NULL), owns_types_(false) {

    validate_dimensions();
    create_datatypes();
    owns_types_ = true;
}

GhostCellExchange::~GhostCellExchange() {
    if (owns_types_) {
        free_datatypes();
    }
}

GhostCellExchange::GhostCellExchange(GhostCellExchange&& other) noexcept
    : nx_(other.nx_), ny_(other.ny_), topology(other.topology),
      row_type_(other.row_type_), col_type_(other.col_type_),
      owns_types_(other.owns_types_) {
    other.owns_types_ = false;
    other.row_type_ = MPI_DATATYPE_NULL;
    other.col_type_ = MPI_DATATYPE_NULL;
}

GhostCellExchange& GhostCellExchange::operator=(GhostCellExchange&& other) noexcept {
    if (this != &other) {
        if (owns_types_) {
            free_datatypes();
        }

        nx_ = other.nx_;
        ny_ = other.ny_;
        // Note: topology is a reference, cannot be reassigned
        row_type_ = other.row_type_;
        col_type_ = other.col_type_;
        owns_types_ = other.owns_types_;

        other.owns_types_ = false;
        other.row_type_ = MPI_DATATYPE_NULL;
        other.col_type_ = MPI_DATATYPE_NULL;
    }
    return *this;
}

void GhostCellExchange::exchange(utils::Array2D& array) {
    if (!validate_array_size(array)) {
        throw std::invalid_argument(
            "Array size does not match expected dimensions. "
            "Expected: " + std::to_string(ny_ + 2) + " x " + std::to_string(nx_ + 2) +
            ", Got: " + std::to_string(array.rows()) + " x " + std::to_string(array.cols())
        );
    }

    const auto& neighbors = topology.neighbors();
    MPI_Comm comm = topology.communicator();
    int my_rank = topology.rank();

    // Array dimensions with ghost cells
    int total_cols = nx_ + 2;
    int total_rows = ny_ + 2;

    // Define row and column indices
    int south_row = 1;              // Interior row adjacent to south ghost
    int north_row = ny_;            // Interior row adjacent to north ghost
    int ghost_south_row = 0;        // South ghost row
    int ghost_north_row = ny_ + 1;   // North ghost row

    int west_col = 1;               // Interior column adjacent to west ghost
    int east_col = nx_;             // Interior column adjacent to east ghost
    int ghost_west_col = 0;         // West ghost column
    int ghost_east_col = nx_ + 1;   // East ghost column

    MPI_Status status;

    // Exchange 1: Send south row to South, receive from North
    // Send: array(south_row, :) to South
    // Receive: array(ghost_north_row, :) from North
    MPI_Sendrecv(
        &array(south_row, 0), 1, row_type_, neighbors.south, 0,
        &array(ghost_north_row, 0), 1, row_type_, neighbors.north, 0,
        comm, &status
    );

    // Exchange 2: Send north row to North, receive from South
    // Send: array(north_row, :) to North
    // Receive: array(ghost_south_row, :) from South
    MPI_Sendrecv(
        &array(north_row, 0), 1, row_type_, neighbors.north, 1,
        &array(ghost_south_row, 0), 1, row_type_, neighbors.south, 1,
        comm, &status
    );

    // Exchange 3: Send west column to West, receive from East
    // Send: array(:, west_col) to West
    // Receive: array(:, ghost_east_col) from East
    MPI_Sendrecv(
        &array(0, west_col), 1, col_type_, neighbors.west, 2,
        &array(0, ghost_east_col), 1, col_type_, neighbors.east, 2,
        comm, &status
    );

    // Exchange 4: Send east column to East, receive from West
    // Send: array(:, east_col) to East
    // Receive: array(:, ghost_west_col) from West
    MPI_Sendrecv(
        &array(0, east_col), 1, col_type_, neighbors.east, 3,
        &array(0, ghost_west_col), 1, col_type_, neighbors.west, 3,
        comm, &status
    );
}

void GhostCellExchange::exchange_async(utils::Array2D& array, MPI_Request* requests) {
    if (!validate_array_size(array)) {
        throw std::invalid_argument(
            "Array size does not match expected dimensions. "
            "Expected: " + std::to_string(ny_ + 2) + " x " + std::to_string(nx_ + 2) +
            ", Got: " + std::to_string(array.rows()) + " x " + std::to_string(array.cols())
        );
    }

    const auto& neighbors = topology.neighbors();
    MPI_Comm comm = topology.communicator();

    // Define row and column indices
    int south_row = 1;
    int north_row = ny_;
    int ghost_south_row = 0;
    int ghost_north_row = ny_ + 1;

    int west_col = 1;
    int east_col = nx_;
    int ghost_west_col = 0;
    int ghost_east_col = nx_ + 1;

    // Post all non-blocking receives first (to allow overlapping)
    // requests[0]: Irecv from North (for south boundary)
    MPI_Irecv(
        &array(ghost_north_row, 0), 1, row_type_, neighbors.north, 0,
        comm, &requests[0]
    );

    // requests[2]: Irecv from South (for north boundary)
    MPI_Irecv(
        &array(ghost_south_row, 0), 1, row_type_, neighbors.south, 1,
        comm, &requests[2]
    );

    // requests[4]: Irecv from East (for west boundary)
    MPI_Irecv(
        &array(0, ghost_east_col), 1, col_type_, neighbors.east, 2,
        comm, &requests[4]
    );

    // requests[6]: Irecv from West (for east boundary)
    MPI_Irecv(
        &array(0, ghost_west_col), 1, col_type_, neighbors.west, 3,
        comm, &requests[6]
    );

    // Post all non-blocking sends
    // requests[1]: Isend to South (from south boundary)
    MPI_Isend(
        &array(south_row, 0), 1, row_type_, neighbors.south, 0,
        comm, &requests[1]
    );

    // requests[3]: Isend to North (from north boundary)
    MPI_Isend(
        &array(north_row, 0), 1, row_type_, neighbors.north, 1,
        comm, &requests[3]
    );

    // requests[5]: Isend to West (from west boundary)
    MPI_Isend(
        &array(0, west_col), 1, col_type_, neighbors.west, 2,
        comm, &requests[5]
    );

    // requests[7]: Isend to East (from east boundary)
    MPI_Isend(
        &array(0, east_col), 1, col_type_, neighbors.east, 3,
        comm, &requests[7]
    );
}

bool GhostCellExchange::validate_array_size(const utils::Array2D& array) const {
    return array.rows() == static_cast<size_t>(ny_ + 2) &&
           array.cols() == static_cast<size_t>(nx_ + 2);
}

void GhostCellExchange::wait_all(MPI_Request* requests) const {
    MPI_Status statuses[8];
    int mpi_result = MPI_Waitall(8, requests, statuses);

    if (mpi_result != MPI_SUCCESS) {
        throw std::runtime_error("MPI_Waitall failed for ghost cell exchange");
    }
}

void GhostCellExchange::create_datatypes() {
    int mpi_result;

    // Create row type: contiguous (nx+2) doubles
    // This represents a full row including ghost cells in that row
    mpi_result = MPI_Type_contiguous(nx_ + 2, MPI_DOUBLE, &row_type_);
    if (mpi_result != MPI_SUCCESS) {
        throw std::runtime_error("Failed to create row MPI data type");
    }
    mpi_result = MPI_Type_commit(&row_type_);
    if (mpi_result != MPI_SUCCESS) {
        throw std::runtime_error("Failed to commit row MPI data type");
    }

    // Create column type: vector of (ny+2) doubles with stride (nx+2)
    // This allows efficient column communication by skipping rows
    mpi_result = MPI_Type_vector(
        ny_ + 2,     // count: number of elements in each block
        1,           // blocklength: number of elements in each block
        nx_ + 2,     // stride: number of elements between start of each block
        MPI_DOUBLE,  // old type
        &col_type_   // new type
    );
    if (mpi_result != MPI_SUCCESS) {
        throw std::runtime_error("Failed to create column MPI data type");
    }
    mpi_result = MPI_Type_commit(&col_type_);
    if (mpi_result != MPI_SUCCESS) {
        throw std::runtime_error("Failed to commit column MPI data type");
    }
}

void GhostCellExchange::free_datatypes() {
    try {
        if (row_type_ != MPI_DATATYPE_NULL) {
            MPI_Type_free(&row_type_);
            row_type_ = MPI_DATATYPE_NULL;
        }
        if (col_type_ != MPI_DATATYPE_NULL) {
            MPI_Type_free(&col_type_);
            col_type_ = MPI_DATATYPE_NULL;
        }
    } catch (...) {
        // Silently ignore errors in destructor
    }
}

void GhostCellExchange::validate_dimensions() const {
    if (nx_ <= 0) {
        throw std::invalid_argument("nx must be positive");
    }
    if (ny_ <= 0) {
        throw std::invalid_argument("ny must be positive");
    }
}
