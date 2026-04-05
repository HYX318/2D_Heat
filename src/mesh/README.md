# Domain Decomposition Module

## Overview

The DomainDecomposition class provides flexible domain decomposition capabilities for parallel computing with MPI. It divides a 2D global domain into subdomains distributed across multiple processes.

## Classes

### DomainDecomposition

Manages the decomposition of a 2D global domain into subdomains for parallel processing.

#### Key Features

- **Automatic decomposition**: Calculates optimal dimensions closest to a square
- **Manual Cartesian decomposition**: User-specified dimensions
- **Non-uniform decomposition**: Allows different subdomain sizes
- **Coordinate transformation**: Maps between global and local indices
- **Boundary detection**: Identifies process boundaries
- **Load balance analysis**: Computes load distribution statistics

#### Constructors

```cpp
// Automatic decomposition
DomainDecomposition(size_t global_nx, size_t global_ny,
                  const CartesianTopology& topology);

// Manual Cartesian decomposition
DomainDecomposition(size_t global_nx, size_t global_ny,
                  size_t dim_x, size_t dim_y,
                  const CartesianTopology& topology);

// Non-uniform decomposition
DomainDecomposition(std::vector<Subdomain> subdomains,
                  const CartesianTopology& topology);
```

#### Main Methods

- `my_subdomain()`: Get current process's subdomain
- `is_in_my_domain(global_i, global_j)`: Check if coordinates are in local domain
- `find_owner_rank(global_i, global_j)`: Find process owning specific coordinates
- `is_on_boundary(Direction dir)`: Check if on boundary
- `compute_load_balance()`: Calculate load balance statistics

#### Usage Example

```cpp
#include "mpi/mpi_context.hpp"
#include "mpi/cartesian_topology.hpp"
#include "mesh/domain_decomposition.hpp"

int main(int argc, char** argv) {
    MPIContext mpi(argc, argv);
    CartesianTopology topology(MPI_COMM_WORLD);

    // Automatic decomposition
    DomainDecomposition decomp(1000, 800, topology);

    // Get local domain
    const auto& my_domain = decomp.my_subdomain();
    std::cout << "My domain: [" << my_domain.start_i << ":"
              << my_domain.start_i + my_domain.nx << "] x ["
              << my_domain.start_j << ":"
              << my_domain.start_j + my_domain.ny << "]" << std::endl;

    // Coordinate transformation
    size_t global_i = my_domain.start_i + 5;
    size_t global_j = my_domain.start_j + 5;
    size_t local_i = my_domain.global_to_local_i(global_i);
    size_t local_j = my_domain.global_to_local_j(global_j);

    // Load balance analysis
    auto stats = decomp.compute_load_balance();
    std::cout << "Imbalance ratio: " << stats.imbalance_ratio << std::endl;

    return 0;
}
```

### Subdomain

Structure containing information about a local subdomain.

#### Fields

- `start_i`: Global starting index in X direction
- `start_j`: Global starting index in Y direction
- `nx`: Local domain size in X
- `ny`: Local domain size in Y
- `global_nx`: Global domain size in X
- `global_ny`: Global domain size in Y

#### Methods

- `global_to_local_i(global_i)`: Convert global X to local X
- `global_to_local_j(global_j)`: Convert global Y to local Y

### LoadBalanceStats

Statistics about load balancing across subdomains.

#### Fields

- `min_load`: Minimum load across all subdomains
- `max_load`: Maximum load across all subdomains
- `avg_load`: Average load across all subdomains
- `imbalance_ratio`: Load imbalance ratio (max - min) / avg

### DecompositionStrategy

Enumeration of decomposition strategies.

- `Automatic`: Automatically compute optimal decomposition
- `Cartesian`: Regular Cartesian grid
- `Manual`: User-specified decomposition

## Implementation Details

### Automatic Dimension Calculation

The `compute_optimal_dimensions()` function finds dimensions that:
1. Minimize the aspect ratio difference (local domain closest to square)
2. Balance load as evenly as possible
3. Are factors of the number of processes

Example:
- 4 processes → [2, 2]
- 6 processes → [2, 3]
- 8 processes → [2, 4]
- 12 processes → [3, 4]

### Cartesian Decomposition

For regular Cartesian decomposition, the global domain is divided into:
- `base_nx = global_nx / dim_x` cells per subdomain in X
- `base_ny = global_ny / dim_y` cells per subdomain in Y
- Remainder cells are distributed among the first `rem_nx` columns and `rem_ny` rows

### Manual Decomposition Validation

The `validate_manual_decomposition()` method ensures:
- All subdomains cover the entire global domain
- No overlapping regions exist
- Every cell is assigned to exactly one subdomain

## Thread Safety

The DomainDecomposition class is not thread-safe for concurrent modification. However, const methods can be called concurrently from multiple threads after construction.

## Exception Safety

All methods are exception-safe and follow RAII principles. The class validates:
- Valid Cartesian topology
- Correct dimension specifications
- Complete and non-overlapping domain coverage

## Dependencies

- `mpi/cartesian_topology.hpp`: For MPI Cartesian topology management
- `<vector>`, `<cmath>`, `<algorithm>`, `<stdexcept>`: Standard library

## Build Instructions

Add to CMakeLists.txt:

```cmake
add_library(heat_equation_mesh STATIC
    mesh/domain_decomposition.cpp
    # ... other mesh sources
)

target_include_directories(heat_equation_mesh
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(heat_equation_mesh
    PUBLIC
        heat_equation_utils
)
```

## See Also

- `cartesian_topology.hpp`: MPI Cartesian topology management
- `mpi_context.hpp`: MPI initialization and management
- Examples in `examples/domain_decomposition_example.cpp`
