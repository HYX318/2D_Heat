# Contributing to 2D Heat Equation Solver

Thank you for your interest in contributing to the 2D Heat Equation Solver project! This document provides guidelines and instructions for contributors.

---

# Table of Contents

1. [Getting Started](#getting-started)
2. [Development Environment Setup](#development-environment-setup)
3. [Code Style Guidelines](#code-style-guidelines)
4. [Testing Requirements](#testing-requirements)
5. [Commit Messages](#commit-messages)
6. [Pull Request Process](#pull-request-process)
7. [Documentation Standards](#documentation-standards)
8. [Issue Reporting](#issue-reporting)

---

# Getting Started

## Prerequisites

Before contributing, ensure you have:

- Git installed and configured
- CMake 3.16 or higher
- C++17 compatible compiler
- MPI (OpenMPI or MPICH)
- Google Test (auto-downloaded)
- A GitHub account

## Forking the Repository

```bash
# Fork the repository on GitHub
# Clone your fork
git clone https://github.com/yourusername/2D_Heat.git
cd 2D_Heat

# Add upstream remote
git remote add upstream https://github.com/originalowner/2D_Heat.git
```

## Creating a Branch

```bash
# Create a new branch for your changes
git checkout -b feature/my-feature

# Or for bug fixes
git checkout -b fix/issue-123

# Or for documentation
git checkout -b docs/update-api-docs
```

---

# Development Environment Setup

## Building in Debug Mode

```bash
# Create build directory
mkdir -p build
cd build

# Configure with debug options
cmake .. \
    -DCMAKE_BUILD_TYPE=Debug \
    -DENABLE_MPI=ON \
    -DBUILD_TESTS=ON \
    -DBUILD_EXAMPLES=ON

# Build
make -j$(nproc)
```

## Building in Release Mode

```bash
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

## IDE Setup

### VS Code

Install recommended extensions:
- C/C++ (Microsoft)
- CMake Tools (Microsoft)
- C/C++ Test Explorer (Microsoft)

Create `.vscode/settings.json`:
```json
{
    "cmake.configureOnOpen": true,
    "C_Cpp.default.configurationProvider": "ms-vscode.cmake-tools",
    "C_Cpp.default.cppStandard": "c++17"
}
```

### CLion

CLion will automatically detect CMakeLists.txt. Configure with:
- Build type: Debug
- C++ standard: C++17
- Compiler: mpicxx

---

# Code Style Guidelines

## C++ Code Style

### Indentation and Formatting

- Use 4 spaces for indentation (no tabs)
- Maximum line length: 120 characters
- One statement per line
- Opening braces on same line as statement (K&R style)

**Example:**
```cpp
// Good
if (condition) {
    do_something();
} else {
    do_something_else();
}

// Bad (tabs)
if (condition)
{
	do_something();
}

// Bad (Allman style)
if (condition)
{
    do_something();
}
else
{
    do_something_else();
}
```

### Naming Conventions

**Classes**: PascalCase
```cpp
class Mesh2D { ... };
class ConjugateGradientSolver { ... };
```

**Functions**: camelCase
```cpp
void compute_laplacian();
double get_solution_value();
```

**Variables**: snake_case
```cpp
double lambda;
size_t num_iterations;
```

**Member Variables**: snake_case with trailing underscore
```cpp
class MyClass {
private:
    double value_;
    size_t count_;
};
```

**Constants**: UPPER_SNAKE_CASE
```cpp
constexpr double PI = 3.14159265358979323846;
constexpr size_t MAX_ITERATIONS = 10000;
```

**Macros**: UPPER_SNAKE_CASE (avoid when possible)
```cpp
#define LOG_ERROR(msg) std::cerr << msg << std::endl
```

### File Organization

- Header files: `.hpp` extension
- Implementation files: `.cpp` extension
- One class per file (header + implementation)
- Organize by module (mesh, solver, mpi, utils)

**Directory Structure:**
```
src/
├── mesh/
│   ├── mesh2d.hpp
│   ├── mesh2d.cpp
│   ├── domain_decomposition.hpp
│   └── domain_decomposition.cpp
├── solver/
│   ├── jacobi_solver.hpp
│   └── jacobi_solver.cpp
└── utils/
    ├── array2d.hpp
    └── logger.hpp/cpp
```

### Comments and Documentation

**File Headers:**
```cpp
/**
 * @file mesh2d.cpp
 * @brief Implementation of 2D mesh with ghost cell support
 *
 * This file implements the Mesh2D class which provides:
 * - 2D rectangular mesh management
 * - Ghost cell support for MPI
 * - Boundary condition handling
 * - Coordinate transformations
 */
```

**Function Documentation:**
```cpp
/**
 * @brief Compute the Laplacian using 5-point stencil
 * @param result Output array for Laplacian
 * @throws std::invalid_argument if result size doesn't match
 *
 * Computes Laplacian using central difference approximation:
 * ∇²u = (u_{i+1,j} - 2u_{i,j} + u_{i-1,j})/hx² +
 *       (u_{i,j+1} - 2u_{i,j} + u_{i,j-1})/hy²
 */
void compute_laplacian(utils::Array2D& result) const;
```

**Class Documentation:**
```cpp
/**
 * @class Mesh2D
 * @brief 2D rectangular mesh with ghost cell support
 *
 * This class manages a 2D rectangular grid with support for:
 * - Ghost cells for MPI domain decomposition
 * - Dirichlet boundary conditions
 * - Coordinate transformations (local/global, grid/physical)
 * - Numerical operations (Laplacian, norms, arithmetic)
 */
```

**Inline Comments:**
```cpp
// Compute Jacobi coefficient
double coeff = 1.0 / (1.0 + 4.0 * params.lambda);

// Update interior points
for (size_t j = 1; j <= ny; ++j) {
    for (size_t i = 1; i <= nx; ++i) {
        // 5-point stencil
        solution(j, i) = coeff * (
            lambda * (solution(j+1, i) + solution(j-1, i) +
                      solution(j, i+1) + solution(j, i-1)) +
            rhs(j, i)
        );
    }
}
```

### C++ Best Practices

**Use RAII:**
```cpp
// Good: Automatic resource management
{
    MPIContext mpi(argc, argv);
    // ... code ...
} // MPI automatically finalized

// Bad: Manual resource management
MPI_Init(&argc, &argv);
// ... code ...
MPI_Finalize(); // Easy to forget
```

**Use smart pointers:**
```cpp
// Good: Automatic memory management
std::unique_ptr<int[]> data(new int[size]);
std::shared_ptr<Solver> solver(new JacobiSolver());

// Bad: Manual memory management
int* data = new int[size];
// ... code ...
delete[] data; // Easy to forget
```

**Use move semantics:**
```cpp
// Good: Efficient transfer
utils::Array2D a(100, 100);
utils::Array2D b = std::move(a);  // Fast, no copy

// Bad: Unnecessary copy
utils::Array2D a(100, 100);
utils::Array2D b = a;  // Slow, copies data
```

**Use constexpr:**
```cpp
// Good: Compile-time constant
constexpr double PI = 3.14159265358979323846;
constexpr size_t MAX_ITER = 10000;

// Bad: Runtime constant
const double PI = 3.14159265358979323846;
const size_t MAX_ITER = 10000;
```

**Use exceptions for error handling:**
```cpp
// Good: Exception-based error handling
try {
    solver.solve(rhs, solution, params);
} catch (const std::runtime_error& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
}

// Bad: Error codes
int error = solver.solve(rhs, solution, params);
if (error != 0) {
    std::cerr << "Error occurred" << std::endl;
    return 1;
}
```

### Modern C++ Features

**Use auto for type deduction:**
```cpp
// Good: Clearer, less verbose
auto it = std::find(vec.begin(), vec.end(), value);
auto result = solver.get_stats();

// Bad: Verbose
std::vector<int>::iterator it = std::find(vec.begin(), vec.end(), value);
SolverStats result = solver.get_stats();
```

**Use range-based for loops:**
```cpp
// Good: Cleaner, less error-prone
for (const auto& elem : container) {
    process(elem);
}

// Bad: Manual index management
for (size_t i = 0; i < container.size(); ++i) {
    process(container[i]);
}
```

**Use lambdas:**
```cpp
// Good: Concise callback
mesh.apply_bc([](double x, double y, double t) {
    return std::sin(x) * std::cos(y);
}, 0.0);

// Bad: Separate function
double bc_func(double x, double y, double t) {
    return std::sin(x) * std::cos(y);
}
mesh.apply_bc(bc_func, 0.0);
```

---

# Testing Requirements

## Unit Tests

All new functionality must include unit tests using Google Test.

### Test File Organization

Create test files in `tests/unit/`:
```cpp
// tests/unit/test_my_new_feature.cpp
#include <gtest/gtest.h>
#include "path/to/my_new_feature.hpp"

TEST(MyNewFeatureTest, Constructor) {
    // Test constructor
}

TEST(MyNewFeatureTest, ComputeValue) {
    // Test value computation
}
```

### Adding Tests to CMake

Add to `tests/unit/CMakeLists.txt`:
```cmake
add_executable(test_my_new_feature
    test_my_new_feature.cpp
)
target_link_libraries(test_my_new_feature
    PRIVATE
        my_new_feature_lib
        gtest_main
)
add_test(NAME test_my_new_feature COMMAND test_my_new_feature)
```

### Test Naming

Use descriptive test names:
```cpp
// Good: Descriptive
TEST(JacobiSolverTest, ConvergesForSimpleSystem)
TEST(Mesh2DTest, BoundaryConditionsAppliedCorrectly)

// Bad: Vague
TEST(JacobiSolverTest, Test1)
TEST(Mesh2DTest, TestBC)
```

### Test Coverage

Aim for at least 80% code coverage:
- Normal execution paths
- Edge cases and boundary conditions
- Error handling and exception cases
- MPI communication (if applicable)

### Running Tests

```bash
# Run all tests
cd build
ctest

# Run tests with coverage
ctest --enable-coverage
```

## Integration Tests

Add integration tests to `tests/integration/`:
```cpp
// tests/integration/test_full_solver.cpp
TEST(IntegrationTest, FullSolverExecution) {
    // Test complete solver workflow
}
```

## Performance Tests

Add performance benchmarks to `tests/performance/`:
```cpp
// tests/performance/benchmark_solvers.cpp
TEST(PerformanceTest, SolverComparison) {
    // Benchmark different solvers
}
```

---

# Commit Messages

## Commit Message Format

Follow [Conventional Commits](https://www.conventionalcommits.org/):

```
<type>(<scope>): <subject>

<body>

<footer>
```

### Types

- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `style`: Code style changes (formatting, etc.)
- `refactor`: Code refactoring
- `perf`: Performance improvements
- `test`: Adding or updating tests
- `chore`: Maintenance tasks
- `revert`: Revert a previous commit

### Examples

**Feature:**
```
feat(solver): Add conjugate gradient solver

- Implement CG algorithm
- Add preconditioning support
- Update documentation
- Add unit tests

Closes #123
```

**Bug Fix:**
```
fix(mpi): Fix memory leak in ghost cell exchange

- Properly free allocated buffers
- Add error handling
- Update tests

Closes #456
```

**Documentation:**
```
docs(api): Update API documentation

- Add missing function documentation
- Fix typos
- Add examples

Related to #789
```

**Performance:**
```
perf(solver): Optimize Jacobi solver loop order

- Reorder loops for better cache locality
- Add vectorization hints
- Benchmark improvement: 15% faster

Closes #101
```

---

# Pull Request Process

## Before Submitting

1. **Update your branch:**
```bash
git fetch upstream
git rebase upstream/master
```

2. **Run tests:**
```bash
cd build
ctest
```

3. **Check code style:**
```bash
# Use clang-format if available
clang-format -i src/**/*.cpp src/**/*.hpp
```

4. **Build in release mode:**
```bash
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

5. **Update documentation:**
- Add documentation for new features
- Update examples if needed
- Update README if public API changed

## Creating Pull Request

1. **Push to the origin:**
```bash
git push origin feature/my-feature
```

2. **Create PR on GitHub:**
- Go to repository on GitHub
- Click "New Pull Request"
- Select your branch
- Fill out PR template

3. **PR Description:**
```markdown
## Summary
Brief description of changes.

## Changes Made
- Change 1
- Change 2
- Change 3

## Testing
- [x] Unit tests pass
- [x] Integration tests pass
- [x] Manual testing completed

## Breaking Changes
List any breaking changes (if none, say "None")

## Checklist
- [x] Code follows style guidelines
- [x] Tests added/updated
- [x] Documentation updated
- [x] No compiler warnings
- [x] Commits follow conventional format
```

## Review Process

1. **Wait for review**: Maintainers will review your PR
2. **Address feedback**: Make requested changes
3. **Update PR**: Push changes to your branch
4. **Approval**: Once approved, PR will be merged

## Merging

Maintainers will squash and merge your commits to keep history clean.

---

## Documentation Standards

## Inline Documentation

Every public function/class must have documentation:
```cpp
/**
 * @brief Brief description
 * @param param1 Description of parameter
 * @return Description of return value
 * @throws std::runtime_error Description of error condition
 */
```

## README Updates

Update README.md for:
- New features
- API changes
- Build system changes
- Performance improvements

## Examples

Add examples for new functionality:
```cpp
// examples/my_new_feature_example.cpp
#include "path/to/my_new_feature.hpp"

int main() {
    // Demonstrate usage
    return 0;
}
```

## API Documentation

Update `docs/api.md` for:
- New classes
- New methods
- Changed interfaces
- Deprecation notices

---

# Issue Reporting

## Bug Reports

When reporting bugs, include:

1. **Environment:**
   - OS version
   - Compiler version
   - CMake version
   - MPI implementation

2. **Steps to reproduce:**
   ```bash
   # Commands to reproduce
   ```

3. **Expected behavior:**
   What should happen

4. **Actual behavior:**
   What actually happens

5. **Error messages:**
   ```plaintext
   Complete error output
   ```

6. **Additional context:**
   - Configuration used
   - Problem size
   - Any relevant logs

## Feature Requests

When requesting features, include:

1. **Use case:**
   What problem does this solve?

2. **Proposed solution:**
   How should it work?

3. **Alternatives considered:**
   Other approaches you've considered

4. **Additional context:**
   Any relevant information

---

# Code of Conduct

## Our Pledge

We are committed to making participation in our project a harassment-free experience for everyone.

## Our Standards

Examples of behavior that contributes to a positive environment:
- Using welcoming and inclusive language
- Being respectful of differing viewpoints
- Gracefully accepting constructive criticism
- Focusing on what is best for the community

Examples of unacceptable behavior:
- The use of sexualized language or imagery
- Personal attacks or insulting comments
- Public or private harassment
- Publishing others' private information

## Our Responsibilities

Maintainers are responsible for clarifying standards and responding to any inappropriate behavior.

## Scope

This code of conduct applies within all project spaces and public spaces.

---

# Getting Help

If you need help contributing:

1. **Read documentation:**
   - [API Documentation](docs/api.md)
   - [Examples](docs/examples.md)
   - [Architecture](docs/architecture.md)

2. **Ask questions:**
   - GitHub Discussions
   - GitHub Issues (using "question" label)

3. **Contact maintainers:**
   - Email: your.email@example.com
   - Slack: #2d-heat-equation

---

# Recognition

Contributors will be recognized in:
- CONTRIBUTORS.md file
- Release notes
- Project website

---

# License

By contributing, you agree that your contributions will be licensed under the MIT License.

---

**Thank you for contributing!**
