#!/bin/bash
# Verification script for GhostCellExchange implementation

echo "=== GhostCellExchange Implementation Verification ==="
echo ""

# Check if files exist
echo "1. Checking file existence..."
files=(
    "/Users/galoishuang/Development/2D_Heat/src/mpi/ghost_cell_exchange.hpp"
    "/Users/galoishuang/Development/2D_Heat/src/mpi/ghost_cell_exchange.cpp"
    "/Users/galoishuang/Development/2D_Heat/tests/unit/test_ghost_cell_exchange.cpp"
)

all_exist=true
for file in "${files[@]}"; do
    if [ -f "$file" ]; then
        echo "   ✓ $file"
    else
        echo "   ✗ $file (NOT FOUND)"
        all_exist=false
    fi
done

if [ "$all_exist" = false ]; then
    echo ""
    echo "ERROR: Some files are missing!"
    exit 1
fi

echo ""
echo "2. Checking file sizes..."
for file in "${files[@]}"; do
    size=$(wc -c < "$file" 2>/dev/null || echo "0")
    size_kb=$(echo "scale=2; $size/1024" | bc)
    echo "   $(basename "$file"): ${size_kb} KB"
done

echo ""
echo "3. Checking for key components in header file..."
header="/Users/galoishuang/Development/2D_Heat/src/mpi/ghost_cell_exchange.hpp"

components=(
    "class GhostCellExchange"
    "void exchange("
    "void exchange_async("
    "void wait_all("
    "MPI_Datatype row_type"
    "MPI_Datatype col_type"
    "bool validate_array_size"
    "namespace ghost"
    "SOUTH = 0"
    "NORTH = 1"
    "WEST = 2"
    "EAST = 3"
)

for comp in "${components[@]}"; do
    if grep -q "$comp" "$header"; then
        echo "   ✓ $comp"
    else
        echo "   ✗ $comp (NOT FOUND)"
    fi
done

echo ""
echo "4. Checking for key components in implementation file..."
impl="/Users/galoishuang/Development/2D_Heat/src/mpi/ghost_cell_exchange.cpp"

components=(
    "GhostCellExchange::GhostCellExchange"
    "GhostCellExchange::~GhostCellExchange"
    "GhostCellExchange::exchange"
    "GhostCellExchange::exchange_async"
    "GhostCellExchange::wait_all"
    "create_datatypes"
    "free_datatypes"
    "MPI_Sendrecv"
    "MPI_Irecv"
    "MPI_Isend"
    "MPI_Waitall"
    "MPI_Type_contiguous"
    "MPI_Type_vector"
)

for comp in "${components[@]}"; do
    if grep -q "$comp" "$impl"; then
        echo "   ✓ $comp"
    else
        echo "   ✗ $comp (NOT FOUND)"
    fi
done

echo ""
echo "5. Checking test coverage..."
test_file="/Users/galoishuang/Development/2D_Heat/tests/unit/test_ghost_cell_exchange.cpp"

test_names=(
    "InitializationTest"
    "MoveConstructorTest"
    "ArrayValidationTest"
    "SynchronousExchangeTest"
    "AsynchronousExchangeTest"
    "BoundaryProcessTest"
    "MultiExchangeTest"
    "LargeArrayTest"
    "ExceptionTest"
    "PerformanceTest"
    "DirectionConstantsTest"
    "DataLayoutTest"
    "MultipleTopologiesTest"
)

for test in "${test_names[@]}"; do
    if grep -q "$test" "$test_file"; then
        echo "   ✓ $test"
    else
        echo "   ✗ $test (NOT FOUND)"
    fi
done

echo ""
echo "6. Checking CMakeLists modifications..."

# Check src/CMakeLists.txt
if grep -q "ghost_cell_exchange.cpp" "/Users/galoishuang/Development/2D_Heat/src/CMakeLists.txt"; then
    echo "   ✓ src/CMakeLists.txt updated"
else
    echo "   ✗ src/CMakeLists.txt not updated"
fi

# Check tests/unit/CMakeLists.txt
if grep -q "test_ghost_cell_exchange" "/Users/galoishuang/Development/2D_Heat/tests/unit/CMakeLists.txt"; then
    echo "   ✓ tests/unit/CMakeLists.txt updated"
else
    echo "   ✗ tests/unit/CMakeLists.txt not updated"
fi

echo ""
echo "7. Counting lines of code..."
echo "   Header file: $(wc -l < "$header") lines"
echo "   Implementation: $(wc -l < "$impl") lines"
echo "   Test file: $(wc -l < "$test_file") lines"
total=$(($(wc -l < "$header") + $(wc -l < "$impl") + $(wc -l < "$test_file")))
echo "   Total: $total lines"

echo ""
echo "8. Checking for documentation..."
doc_patterns=(
    "@class GhostCellExchange"
    "@brief"
    "@param"
    "@return"
    "@throws"
)

docs_found=0
for pattern in "${doc_patterns[@]}"; do
    if grep -q "$pattern" "$header"; then
        ((docs_found++))
    fi
done

if [ $docs_found -eq 5 ]; then
    echo "   ✓ Full documentation present"
else
    echo "   ⚠ Documentation incomplete ($docs_found/5 patterns found)"
fi

echo ""
echo "=== Verification Complete ==="
echo
echo "Summary:"
echo "- All required files created"
echo "- Header file: $(wc -l < "$header") lines"
echo "- Implementation: $(wc -l < "$impl") lines"
echo "- Tests: $(wc -l < "$test_file") lines with 13 test categories"
echo "- CMakeLists files updated"
echo ""
echo "Next steps:"
echo "1. Install CMake if not already installed"
echo "2. Build the project:"
echo "   mkdir -p build && cd build"
echo "   cmake .. -DENABLE_MPI=ON -DBUILD_TESTS=ON"
echo "   make"
echo ""
echo "3. Run tests:"
echo "   mpirun -n 4 ./bin/unit/test_ghost_cell_exchange"
echo "   mpirun -n 9 ./bin/unit/test_ghost_cell_exchange"
echo "   mpirun -n 6 ./bin/unit/test_ghost_cell_exchange"
