# Legacy Code Compatibility Layer

This directory contains the legacy compatibility layer for the 2D Heat Equation solver.

## Purpose

The compatibility layer allows existing code that uses the old API to continue working while transitioning to the new architecture. All legacy functions are marked with `[[deprecated]]` to encourage migration.

## Components

### Files

- `heat_utils_compat.hpp/cpp` - Compatibility wrappers for HeatUtils functions
- `interfaces_compat.hpp/cpp` - Compatibility wrappers for MPI communication functions
- `legacy_main.cpp` - Example main program using legacy API
- `test_legacy_compat.cpp` - Tests for compatibility layer
- `CMakeLists.txt` - Build configuration
- `MIGRATION_GUIDE.md` - Guide for migrating to new API

## Building

To build the compatibility layer:

```bash
cd build
cmake ..
make heat_legacy heat_equation_legacy_compat
```

## Running

To run the legacy main program:

```bash
# Requires 9 MPI processes for 3x3 decomposition
mpirun -np 9 ./bin/heat_legacy
```

To run the compatibility tests:

```bash
# Serial test
mpirun -np 1 ./bin/test_legacy_compat

# Parallel test (requires 4 processes)
mpirun -np 4 ./bin/test_legacy_compat
```

## Deprecated Functions

The following functions are deprecated and should be replaced:

### HeatUtils Functions

- `ReadParam()` - Use new configuration system
- `Contiguous2D()` - Use `utils::Array2D` or `Mesh2D`
- `Init()` - Use `Mesh2D::apply_bc()`
- `ExactSol()` - Use analytic solution function with `Mesh2D`
- `Analytic()` - Use new analytic solution function
- `Error()` - Use `Array2D` arithmetic
- `TwoNorm()` - Use `Array2D::l2_norm()`
- `InftyNorm()` - Use `Array2D::linfty_norm()`
- `Copy()` - Use `Array2D::copy_from()`
- `Jacobi()` - Use `JacobiSolver<false>`
- `Export()` - Use new I/O utilities

### Interfaces Functions

- `MPIJacobi()` - Use `JacobiSolver<true>`
- `Interfaces()` - Use `GhostCellExchange::exchange()`

## Migration

See [MIGRATION_GUIDE.md](MIGRATION_GUIDE.md) for detailed instructions on migrating to the new API.

## Warnings

When compiling code that uses the legacy compatibility layer, you will see deprecation warnings:

```
warning: 'ReadParam' is deprecated: Use new configuration system. See MIGRATION_GUIDE.md [-Wdeprecated-declarations]
```

These warnings are intentional to encourage migration to the new API.

## Performance

The compatibility layer is designed to maintain the same performance as the original legacy code while using the new architecture internally. However, for best performance, migrate to the new API directly.

## Support

The compatibility layer will be maintained until version 2.0.0, at which point it will be removed. All users should migrate to the new API before then.
