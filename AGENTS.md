# Repository Guidelines

## Project Structure & Module Organization

This is a C++17/CMake project for a 2D heat equation solver with MPI support.

- `src/` contains production code.
- `src/core/solver/` contains Jacobi, SOR, and Conjugate Gradient solvers.
- `src/core/scheme/` contains time integration schemes.
- `src/mesh/` contains mesh, coordinate, and domain decomposition code.
- `src/mpi/` contains MPI context, topology, ghost cell exchange, reductions, and profiling.
- `src/io/` contains parameter parsing and solution export.
- `tests/unit/` contains Google Test unit tests; `tests/performance/` contains benchmarks.
- `examples/` contains runnable examples.
- `visualization/` contains Python tools for running the solver and plotting outputs.
- `docs/` contains API, architecture, examples, and performance documentation.

## Build, Test, and Development Commands

Run commands from the repository root unless noted.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTS=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

Use Release builds for performance checks:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=ON
cmake --build build
```

Run MPI-aware tests directly:

```bash
mpirun -n 4 ./build/bin/heat_equation_unit_tests
```

Run and visualize a simulation:

```bash
cd visualization
python3 run_and_visualize.py --preset small --save-only
```

## Coding Style & Naming Conventions

Use 4-space indentation, C++17, RAII, and exception-safe resource ownership. Keep headers focused and prefer existing module patterns over new abstractions. Classes use `PascalCase` such as `Mesh2D`; functions and variables generally use `snake_case`; private members use a trailing underscore, for example `omega_`.

## Testing Guidelines

Tests use Google Test via CMake. Add or update tests for solver behavior, mesh indexing, MPI guard paths, and numerical residual semantics when changing those areas. Name test files as `tests/unit/test_<component>.cpp`. Normal CTest skips MPI-only paths unless MPI is initialized; use `mpirun` for those.

## Commit & Pull Request Guidelines

Recent history uses short conventional-style prefixes such as `fix:`, `docs:`, and `test:`. Prefer concise messages like `fix: correct SOR relaxation coefficient`. Pull requests should include a clear summary, affected modules, test commands run, and notes for MPI or visualization behavior. Link related issues when available.

## Security & Configuration Tips

Do not commit generated solver outputs, build directories, or visualization artifacts. Keep output paths inside the repository or an approved temporary directory, and avoid directory traversal in visualization configuration.
