#!/bin/bash

# Verification script for CMake build system
# This script checks if all required CMake files are present and well-formed

echo "======================================"
echo "CMake Build System Verification"
echo "======================================"
echo

# Check if CMake is installed
if ! command -v cmake &> /dev/null; then
    echo "ERROR: CMake is not installed or not in PATH"
    echo "Please install CMake:"
    echo "  macOS: brew install cmake"
    echo "  Ubuntu/Debian: sudo apt-get install cmake"
    echo "  Fedora/RHEL: sudo dnf install cmake"
    exit 1
fi

echo "CMake version: $(cmake --version | head -n 1)"
echo

# Check for all required CMake files
echo "Checking CMake files..."
required_files=(
    "CMakeLists.txt"
    "cmake/config.hpp.in"
    "src/CMakeLists.txt"
    "src/utils/CMakeLists.txt"
    "src/mpi/CMakeLists.txt"
    "tests/CMakeLists.txt"
    "tests/unit/CMakeLists.txt"
)

all_files_exist=true
for file in "${required_files[@]}"; do
    if [ -f "$file" ]; then
        echo "  ✓ $file"
    else
        echo "  ✗ $file (MISSING)"
        all_files_exist=false
    fi
done

echo

if [ "$all_files_exist" = false ]; then
    echo "ERROR: Some required CMake files are missing"
    exit 1
fi

# Check MPI availability
echo "Checking MPI..."
if command -v mpicxx &> /dev/null; then
    echo "  ✓ MPI compiler found: $(mpicxx --version | head -n 1)"
elif command -v mpicc &> /dev/null; then
    echo "  ✓ MPI compiler found: $(mpicc --version | head -n 1)"
else
    echo "  ⚠ MPI compiler not found (mpicxx or mpicc)"
    echo "    To install MPI:"
    echo "      macOS: brew install open-mpi"
    echo "      Ubuntu/Debian: sudo apt-get install libopenmpi-dev"
    echo "      Fedora/RHEL: sudo dnf install openmpi-devel"
fi

echo

# Try to configure CMake
echo "Attempting CMake configuration..."
build_dir="build_verify"
mkdir -p "$build_dir"

if cmake -S . -B "$build_dir" 2>&1 | tee cmake_config.log; then
    echo
    echo "  ✓ CMake configuration successful"
    echo
    echo "Build targets (from CMake):"
    cmake --build "$build_dir" --target help 2>&1 | grep -E "^\.\.\. [a-z_]+" | head -20
    echo

    # Clean up
    rm -rf "$build_dir"
    rm cmake_config.log

    echo "======================================"
    echo "✓ Verification Complete!"
    echo "======================================"
    echo
    echo "To build the project:"
    echo "  mkdir build && cd build"
    echo "  cmake .."
    echo "  make -j\$(nproc)"
    echo
    echo "To run tests:"
    echo "  cd build && ctest"
    exit 0
else
    echo
    echo "  ✗ CMake configuration failed"
    echo
    echo "Check cmake_config.log for details"
    echo
    # Don't clean up on failure so user can inspect
    exit 1
fi
