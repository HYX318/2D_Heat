#!/bin/bash
# Simple compilation test for GhostCellExchange

set -e

# Find MPI compiler
if [ -z "$MPICC" ]; then
    MPICC=$(which mpicc 2>/dev/null || echo "")
fi

if [ -z "$MPICC" ]; then
    echo "Error: mpicc not found. Please install MPI."
    exit 1
fi

echo "Using MPI compiler: $MPICC"

# Find C++ compiler (prefer mpic++ for MPI)
CXX=$(which mpic++ 2>/dev/null || which g++ 2>/dev/null || which clang++ 2>/dev/null || echo "")

if [ -z "$CXX" ]; then
    echo "Error: C++ compiler not found."
    exit 1
fi

echo "Using C++ compiler: $CXX"

# Compile source files
echo "Compiling GhostCellExchange source files..."

$CXX -c -std=c++17 -I. -Isrc src/mpi/ghost_cell_exchange.cpp -o /tmp/ghost_cell_exchange.o 2>&1 || {
    echo "Compilation failed for ghost_cell_exchange.cpp"
    exit 1
}

echo "Success: ghost_cell_exchange.cpp compiled successfully"
echo ""
echo "Note: Full compilation requires CMake with Google Test."
echo "To compile and run tests, install CMake and run:"
echo "  mkdir -p build && cd build"
echo "  cmake .."
echo "  make"
echo "  ctest --output-on-failure"
echo ""
echo "To run GhostCellExchange tests with MPI:"
echo "  mpirun -n 4 ./bin/unit/test_ghost_cell_exchange"
