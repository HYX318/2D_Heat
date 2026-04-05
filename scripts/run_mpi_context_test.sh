#!/bin/bash
# Run script for MPIContext unit tests

set -e

NUM_PROCES=4

if [ ! -z "$1" ]; then
    NUM_PROCES=$1
fi

echo "Running MPIContext unit tests with $NUM_PROCES processes..."
echo ""

# Set TMPDIR to /tmp to avoid shared memory errors on macOS
export TMPDIR=/tmp

# Run the tests
mpirun -np $NUM_PROCES ./build/test_mpi_context

echo ""
echo "Tests completed!"
