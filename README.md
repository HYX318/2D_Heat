# 2D Heat Equation Solver

A high-performance parallel computing application for solving the 2D heat diffusion equation using MPI (Message Passing Interface).

---

# Features

- **Modern C++17**: Uses RAII, move semantics, and exception safety
- **Parallel Computing**: MPI-based domain decomposition for multi-node execution
- **Multiple Solvers**: Jacobi, SOR, Conjugate Gradient, and Preconditioned CG
- **Time Integration**: Implicit Euler and Crank-Nicolson schemes
- **Flexible Mesh**: Uniform and non-uniform grids with ghost cell support
- **Performance Optimized**: Communication/computation overlap, vectorization, and cache optimization
- **Comprehensive Testing**: Unit tests with Google Test framework
- **Clean API**: Well-documented interfaces with examples

---

# Quick Start

## Prerequisites

- CMake 3.16 or higher
- C++17 compatible compiler (GCC 7+, Clang 5+, Apple Clang 10+)
- MPI (OpenMPI or MPICH)
- Google Test (optional, auto-downloaded)

## Building

```bash
# Clone repository
git clone https://github.com/yourusername/2D_Heat.git
cd 2D_Heat

# Configure and build
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

# Run tests
ctest
```

## Running Examples

```bash
# Serial example
./build/bin/heat_solver

# Parallel example with 4 processes
mpirun -n 4 ./build/bin/heat_solver_parallel
```

---

# Installation

## System Installation

```bash
cd build
sudo make install
```

This installs:
- Libraries to `/usr/local/lib`
- Headers to `/usr/local/include/2D_Heat`
- Executables to `/usr/local/bin`

## Custom Installation

```bash
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/opt/2D_Heat
sudo make install
```

---

# Usage

## Basic Example

```cpp
#include "mesh/mesh2d.hpp"
#include "core/solver/jacobi_solver.hpp"
#include "mpi/mpi_context.hpp"

int main(int argc, char** argv) {
    // Initialize MPI
    MPIContext mpi(argc, argv);

    // Create mesh (100×100 grid on unit square)
    Mesh2D mesh(100, 100, 1.0, 1.0);

    // Apply boundary conditions
    mesh.apply_dirichlet_bc(0.0);

    // Set initial condition
    mesh(50, 50) = 1.0;

    // Create solver
    JacobiSolver<false> solver;
    SolverParams params(1e-6, 10000, 10, true, 0.25);

    // Solve
    utils::Array2D rhs = mesh.data();
    utils::Array2D solution(100, 100);
    solver.solve(rhs, solution, params);

    return 0;
}
```

## Parallel Example

```cpp
#include "mpi/mpi_context.hpp"
#include "mpi/cartesian_topology.hpp"
#include "mpi/ghost_cell_exchange.hpp"
#include "mesh/mesh2d.hpp"
#include "core/solver/conjugate_gradient_solver.hpp"

int main(int argc, char** argv) {
    // Initialize MPI
    MPIContext mpi(argc, argv);

    // Create Cartesian topology
    CartesianTopology topology(MPI_COMM_WORLD);

    // Create mesh with ghost cells
    Mesh2D mesh(200, 200, 1.0, 1.0, topology);

    // Apply boundary conditions
    mesh.apply_dirichlet_bc(0.0);

    // Create ghost cell exchange
    GhostCellExchange exchange(mesh.nx(), mesh.ny(), topology);

    // Create parallel CG solver
    ConjugateGradientSolver solver(true, 0.25);
    SolverParams params(1e-8, 1000);
    solver.set_restart_threshold(100);

    // Solve
    utils::Array2D rhs = mesh.data();
    utils::Array2D solution = mesh.data();
    solver.solve(rhs, solution, params);

    // Print results on root
    if (mpi.is_root()) {
        auto stats = solver.get_stats();
        std::cout << "Iterations: " << stats.iterations << std::endl;
        std::cout << "Residual: " << stats.residual << std::endl;
    }

    return 0;
}
```

---

# Documentation

## Available Documentation

- **[Architecture Documentation](docs/architecture.md)** - System design and architecture
- **[API Documentation](docs/api.md)** - Complete API reference
- **[Performance Guide](docs/performance_guide.md)** - Optimization and tuning
- **[Examples Guide](docs/examples.md)** - Complete code examples
- **[Migration Guide](docs/migration_guide.md)** - Migrating from legacy code
- **[CMake Guide](README_CMAKE.md)** - Build system details
- **[Quick Start](QUICKSTART.md)** - Getting started guide

---

# Project Structure

```
2D_Heat/
├── src/                          # Source code
│   ├── core/                      # Core solver and time stepping
│   │   ├── solver/                # Iterative solvers
│   │   │   ├── jacobi_solver.hpp/cpp
│   │   │   ├── sor_solver.hpp/cpp
│   │   │   └── conjugate_gradient_solver.hpp/cpp
│   │   └── timestepping/         # Time integration schemes
│   │       ├── implicit_euler.hpp/cpp
│   │       └── crank_nicolson.hpp/cpp
│   ├── mesh/                      # Mesh and domain management
│   │   ├── mesh2d.hpp/cpp
│   │   ├── domain_decomposition.hpp/cpp
│   │   └── coordinate_system.hpp/cpp
│   ├── mpi/                       # MPI infrastructure
│   │   ├── mpi_context.hpp/cpp
│   │   ├── cartesian_topology.hpp/cpp
│   │   ├── ghost_cell_exchange.hpp/cpp
│   │   ├── reduction_ops.hpp/cpp
│   │   └── profiler.hpp/cpp
│   └── utils/                     # Utility classes
│       ├── array2d.hpp
│       ├── logger.hpp/cpp
│       ├── timer.hpp
│       ├── param_reader.hpp.hpp
│       └── solution_exporter.hpp.cpp
├── tests/                         # Tests
│   ├── unit/                      # Unit tests
│   └── performance/               # Performance benchmarks
├── examples/                      # Example programs
├── docs/                          # Documentation
│   ├── architecture.md
│   ├── api.md
│   ├── performance_guide.md
│   ├── examples.md
│   └── migration_guide.md
├── CMakeLists.txt                # Main CMake configuration
├── README.md                      # This file
└── QUICKSTART.md                  # Quick start guide
```

---

# Build Options

## CMake Options

```bash
# Build type
cmake .. -DCMAKE_BUILD_TYPE=Release    # Optimized build
cmake .. -DCMAKE_BUILD_TYPE=Debug      # Debug build

# MPI support
cmake .. -DENABLE_MPI=ON              # Enable MPI (default)
cmake .. -DENABLE_MPI=OFF             # Disable MPI

# Build targets
cmake .. -DBUILD_TESTS=ON              # Build tests (default)
cmake .. -DBUILD_EXAMPLES=ON          # Build examples (default)

# Compiler
cmake .. -DCMAKE_CXX_COMPILER=mpicxx   # Use MPI compiler
```

## Make Targets

```bash
# Build all
make -j$(nproc)

# Build specific target
make heat_solver
make test_jacobi_solver

# Install
sudo make install

# Clean
make clean
make distclean
```

---

# Testing

## Run All Tests

```bash
cd build
ctest
```

## Run Tests Verbosely

```bash
ctest --verbose
```

## Run Specific Test

```bash
ctest -R Array2D
ctest -R Jacobi
ctest -R MPI
```

## Run with MPI

```bash
mpirun -n 4 ./bin/test_ghost_cell_exchange
mpirun -n 8 ./bin/test_conjugate_gradient_solver
```

---

# Performance

## Solver Performance Comparison

| Solver | Grid Size | Iterations | Time (s) | Speedup |
|--------|------------|------------|----------|----------|
| Jacobi | 100×100 | 18,500 | 2.45 | 1.0× |
| SOR | 100×100 | 850 | 0.32 | 7.7× |
| CG | 100×100 | 180 | 0.15 | 16.3× |
| PCG | 100×100 | 120 | 0.12 | 20.4× |

## Strong Scaling (500×500 Grid)

| Processes | Time (s) | Speedup | Efficiency |
|-----------|----------|----------|-------------|
| 1 | 145.2 | 1.00 | 100% |
| 2 | 76.8 | 1.89 | 94.5% |
| 4 | 41.2 | 3.53 | 88.3% |
| 8 | 23.5 | 6.18 | 77.3% |
| 16 | 13.8 | 10.52 | 65.8% |

## Optimization Tips

1. **Use Appropriate Solver**: CG/PCG for large problems, SOR for medium problems
2. **Enable Optimizations**: Build with `-DCMAKE_BUILD_TYPE=Release`
3. **Process Count**: Use one process per CPU core
4. **Communication Overlap**: Automatic with GhostCellExchange
5. **Preconditioning**: Enable for ill-conditioned systems

---

# System Requirements

## Minimum Requirements

- **Compiler**: GCC 7+ or Clang 5+
- **CMake**: 3.16 or higher
- **MPI**: OpenMPI 2.0+ or MPICH 3.0+
- **RAM**: 1 GB minimum (depends on problem size)
- **Disk**: 100 MB for installation

## Recommended Configuration

- **Compiler**: GCC 9+ or Clang 10+
- **CMake**: 3.20 or higher
- **MPI**: OpenMPI 4.0+
- **RAM**: 4 GB or more
- **CPU**: Multi-core processor for parallel execution
- **Disk**: 1 GB for installation and examples

## OS Support

- **Linux**: Ubuntu 18.04+, Fedora 30+, CentOS 7+
- **macOS**: 10.14+ (Mojave or later)
- **Windows**: WSL (Windows Subsystem for Linux)

---

# Troubleshooting

## MPI Not Found

```bash
# Specify MPI compiler explicitly
cmake .. -DCMAKE_CXX_COMPILER=mpicxx

# Or find MPI manually
cmake .. -DMPI_CXX_COMPILER=/usr/lib/openmpi/bin/mpicxx
```

## Build Errors

```bash
# Clean and rebuild
rm -rf build
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j$(nproc)
```

## Test Failures

```bash
# Run tests with verbose output
ctest --verbose

# Check specific test
ctest -R TestName --output-on-failure

# Run with valgrind for memory errors
ctest -T memcheck
```

## Performance Issues

1. Check process count (one per core)
2. Verify build type (Release mode)
3. Enable compiler optimizations
4. Profile with built-in Profiler
5. Consider different solver

---

# Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## Development Workflow

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Run tests and ensure they pass
6. Update documentation
7. Submit a pull request

## Code Style

- C++17 standard
- RAII pattern for resource management
- Exception safety
- Clear, self-documenting code
- Comprehensive comments
- Follow existing code style

---

# License

This project is licensed under the MIT License - see LICENSE file for details.

---

# Citation

If you use this software in your research, please cite:

```
@software{2d_heat_equation_solver,
  author = {Your Name},
  title = {2D Heat Equation Solver},
  year = {2024},
  url = {https://github.com/yourusername/2D_Heat}
}
```

---

# Acknowledgments

- MPI Forum for MPI standard
- Google Test for testing framework
- CMake for build system
- The C++ community for best practices

---

# Contact

- **Issues**: [GitHub Issues](https://github.com/yourusername/2D_Heat/issues)
- **Discussions**: [GitHub Discussions](https://github.com/yourusername/2D_Heat/discussions)
- **Email**: your.email@example.com

---

# Links

- [Documentation](docs/)
- [Examples](examples/)
- [API Reference](docs/api.md)
- [Performance Guide](docs/performance_guide.md)
- [Architecture](docs/architecture.md)
- [Migration Guide](docs/migration_guide.md)

---

**Happy Computing!**
