#!/bin/bash

###############################################################################
# run_all_tests.sh - Comprehensive test runner for 2D Heat Equation Solver
#
# This script automates the entire testing process:
#   - Build all tests
#   - Run non-MPI tests
#   - Run MPI tests with specified number of processes
#   - Generate test reports
#
# Usage:
#   ./run_all_tests.sh [options]
#
# Options:
#   -b, --build-only      Only build tests, don't run them
#   -r, --run-only        Only run tests (skip build)
#   -n, --num-procs N     Number of MPI processes (default: 4)
#   -v, --verbose         Enable verbose output
#   -c, --clean           Clean build artifacts before building
#   -h, --help            Show this help message
###############################################################################

set -e  # Exit on error (but will be overridden for test failures)

# Default values
NUM_PROCS=4
BUILD_ONLY=false
RUN_ONLY=false
VERBOSE=false
CLEAN=false
BUILD_DIR="build"
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Print colored message
print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Print header
print_header() {
    echo ""
    echo "=========================================="
    echo "  2D Heat Equation Solver Test Suite"
    echo "=========================================="
    echo "Project Root: $PROJECT_ROOT"
    echo "Build Directory: $BUILD_DIR"
    echo "Number of MPI Processes: $NUM_PROCS"
    echo "=========================================="
    echo ""
}

# Parse command line arguments
parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -b|--build-only)
                BUILD_ONLY=true
                shift
                ;;
            -r|--run-only)
                RUN_ONLY=true
                shift
                ;;
            -n|--num-procs)
                NUM_PROCS="$2"
                shift 2
                ;;
            -v|--verbose)
                VERBOSE=true
                set -x  # Enable command tracing
                shift
                ;;
            -c|--clean)
                CLEAN=true
                shift
                ;;
            -h|--help)
                show_help
                exit 0
                ;;
            *)
                print_error "Unknown option: $1"
                show_help
                exit 1
                ;;
        esac
    done
}

# Show help message
show_help() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  -b, --build-only      Only build tests, don't run them"
    echo "  -r, --run-only        Only run tests (skip build)"
    echo "  -n, --num-procs N     Number of MPI processes (default: 4)"
    echo "  -v, --verbose         Enable verbose output"
    echo "  -c, --clean           Clean build artifacts before building"
    echo "  -h, --help            Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0                    # Build and run all tests"
    echo "  $0 -n 8               # Run with 8 MPI processes"
    echo "  $0 -b                 # Only build tests"
    echo "  $0 -r                 # Only run tests (skip build)"
    echo "  $0 -c -n 4            # Clean, build, and run with 4 processes"
}

# Check prerequisites
check_prerequisites() {
    print_info "Checking prerequisites..."

    # Check if CMake is installed
    if ! command -v cmake &> /dev/null; then
        print_error "CMake is not installed. Please install CMake."
        exit 1
    fi
    print_success "CMake found: $(cmake --version | head -n 1)"

    # Check if MPI is installed
    if ! command -v mpic++ &> /dev/null; then
        print_error "MPI is not installed. Please install OpenMPI or MPICH."
        exit 1
    fi

    # Check if mpirun is available
    if ! command -v mpirun &> /dev/null && ! command -v mpiexec &> /dev/null; then
        print_error "Neither mpirun nor mpiexec found. Cannot run MPI tests."
        exit 1
    fi

    # Detect MPI launcher
    if command -v mpirun &> /dev/null; then
        MPI_LAUNCHER="mpirun"
        MPI_VERSION=$(mpirun --version 2>&1 | head -n 1)
    else
        MPI_LAUNCHER="mpiexec"
        MPI_VERSION=$(mpiexec --version 2>&1 | head -n 1)
    fi
    print_success "MPI launcher: $MPI_LAUNCHER ($MPI_VERSION)"

    print_success "All prerequisites satisfied"
}

# Clean build directory
clean_build() {
    if [ "$CLEAN" = true ]; then
        print_info "Cleaning build directory..."
        rm -rf "$BUILD_DIR"
        print_success "Build directory cleaned"
    fi
}

# Create build directory and configure with CMake
configure_cmake() {
    print_info "Configuring project with CMake..."

    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"

    # Run CMake configuration
    if ! cmake .. -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTS=ON; then
        print_error "CMake configuration failed"
        exit 1
    fi

    print_success "CMake configuration completed"
}

# Build all tests
build_tests() {
    print_info "Building all tests..."

    cd "$BUILD_DIR"

    # Build tests
    if ! cmake --build . --target build_tests; then
        print_error "Build failed"
        exit 1
    fi

    print_success "All tests built successfully"
}

# Run non-MPI tests
run_non_mpi_tests() {
    print_info "Running non-MPI tests..."

    cd "$BUILD_DIR"

    # Run Array2D tests
    print_info "Running Array2D tests..."
    if ./bin/test_array2d; then
        print_success "Array2D tests passed"
    else
        print_error "Array2D tests failed"
        return 1
    fi
}

# Run MPI tests
run_mpi_tests() {
    print_info "Running MPI tests with $NUM_PROCS processes..."

    cd "$BUILD_DIR"

    # Set TMPDIR to avoid shared memory errors on macOS
    export TMPDIR=/tmp

    # Run MPI context tests
    print_info "Running MPIContext tests..."
    if $MPI_LAUNCHER -n $NUM_PROCS ./bin/test_mpi_context; then
        print_success "MPIContext tests passed"
    else
        print_error "MPIContext tests failed"
        return 1
    fi
}

# Run CTest (if available)
run_ctest() {
    if command -vctest &> /dev/null; then
        print_info "Running CTest..."

        cd "$BUILD_DIR"

        # Run non-MPI tests with CTest
        print_info "Running non-MPI tests via CTest..."
        ctest --output-on-failure -R "Array2D" || {
            print_warning "Some CTest tests failed"
        }

        # Run MPI tests with CTest
        print_info "Running MPI tests via CTest..."
        ctest --output-on-failure -R "MPI" || {
            print_warning "Some CTest tests failed"
        }
    else
        print_warning "CTest not available, skipping"
    fi
}

# Generate test report
generate_report() {
    print_info "Generating test report..."

    REPORT_FILE="$BUILD_DIR/test_report.txt"

    cat > "$REPORT_FILE" << EOF
==========================================
  2D Heat Equation Solver Test Report
==========================================
Generated: $(date)
Project Root: $PROJECT_ROOT
Build Directory: $BUILD_DIR
Number of MPI Processes: $NUM_PROCS
==========================================

Environment:
- CMake Version: $(cmake --version | head -n 1)
- C++ Compiler: $(cmake --build . --verbose 2>&1 | grep "CXX_COMPILER" | cut -d'=' -f2 | head -n 1)
- MPI Launcher: $MPI_LAUNCHER

Test Executables:
- test_array2d: $([ -f "$BUILD_DIR/bin/test_array2d" ] && echo "Available" || echo "Not found")
- test_mpi_context: $([ -f "$BUILD_DIR/bin/test_mpi_context" ] && echo "Available" || echo "Not found")

==========================================
End of Report
==========================================
EOF

    print_success "Test report generated: $REPORT_FILE"
    cat "$REPORT_FILE"
}

# Main function
main() {
    parse_args "$@"
    print_header
    check_prerequisites

    if [ "$RUN_ONLY" = false ]; then
        clean_build
        configure_cmake
        build_tests
    fi

    if [ "$BUILD_ONLY" = false ]; then
        # Track overall test result
        TEST_RESULT=0

        # Run non-MPI tests
        if ! run_non_mpi_tests; then
            TEST_RESULT=1
        fi

        # Run MPI tests
        if ! run_mpi_tests; then
            TEST_RESULT=1
        fi

        # Run CTest (optional)
        run_ctest

        # Generate report
        generate_report

        # Print final result
        echo ""
        if [ $TEST_RESULT -eq 0 ]; then
            print_success "=========================================="
            print_success "  ALL TESTS PASSED!"
            print_success "=========================================="
        else
            print_error "=========================================="
            print_error "  SOME TESTS FAILED"
            print_error "=========================================="
            exit 1
        fi
    fi
}

# Run main function
main "$@"
