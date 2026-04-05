#!/bin/bash
# Build script for MPIContext unit tests

set -e

echo "Building MPIContext unit tests..."

# Create build directory if it doesn't exist
"mkdir" -p build

# Compile MPIContext
echo "Compiling MPIContext..."
mpic++ -Wall -g -c src/mpi/mpi_context.cpp -o build/mpi_context.o

# Compile the test
echo "Compiling test_mpi_context..."
mpic++ -Wall -g -std=c++11 build/mpi_context.o tests/unit/test_mpi_context.cpp -o build/test_mpi_context -lgtest -lgtest_main -lpthread

echo "Build complete!"
echo "Executable: build/test_mpi_context"
echo ""
echo "To run tests with 4 processes:"
echo "  mpirun -np 4 ./build/test_mpi_context"
