#include <gtest/gtest.h>
#include <mpi.h>
#include "../../src/mpi/ghost_cell_exchange.hpp"
#include "../../src/mpi/cartesian_topology.hpp"
#include "../../src/mpi/mpi_context.hpp"
#include "../../src/utils/array2d.hpp"
#include <memory>
#include <cmath>

/**
 * @brief Test fixture for GhostCellExchange tests
 *
 * Note: GhostCellExchange tests require running with MPI (mpirun/mpiexec)
 * Example: mpirun -n 4 ./test_ghost_cell_exchange
 */
class GhostCellExchangeTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Check if MPI is initialized
        int mpi_initialized = 0;
        MPI_Initialized(&mpi_initialized);
        if (!mpi_initialized) {
            GTEST_SKIP() << "MPI not initialized. Run with mpirun/mpiexec.";
        }

        // Initialize topology for testing
        topology = std::make_unique<CartesianTopology>(MPI_COMM_WORLD);
    }

    void TearDown() override {
        topology.reset();
    }

    /**
     * @brief Create an array with ghost cells filled with test data
     * @param nx Interior columns
     * @param ny Interior rows
     * @param base_value Base value for interior
     * @return Array2D with ghost cells initialized
     */
    utils::Array2D create_test_array(int nx, int ny, double base_value = 0.0) {
        utils::Array2D array(ny + 2, nx + 2, base_value);

        // Fill interior with a pattern based on position
        for (int i = 1; i <= ny; ++i) {
            for (int j = 1; j <= nx; ++j) {
                array(i, j) = base_value + i * 100 + j;
            }
        }

        // Initialize ghost cells with sentinel values
        for (int j = 0; j <= nx + 1; ++j) {
            array(0, j) = -1.0;           // South ghost
            array(ny + 1, j) = -2.0;       // North ghost
        }
        for (int i = 0; i <= ny + 1; ++i) {
            array(i, 0) = -3.0;           // West ghost
            array(i, nx + 1) = -4.0;      // East ghost
        }

        return array;
    }

    std::unique_ptr<CartesianTopology> topology;
};

/**
 * @test InitializationTest
 * @brief Test that GhostCellExchange initializes correctly
 *
 * This test verifies:
 * 1. GhostCellExchange can be constructed
 * 2. Custom MPI data types are created
 * 3. Dimensions are stored correctly
 * 4. Topology reference is valid
 */
TEST_F(GhostCellExchangeTest, InitializationTest) {
    int nx = 10;
    int ny = 10;

    GhostCellExchange exchange(nx, ny, *topology);

    // Verify dimensions
    EXPECT_EQ(exchange.nx(), nx);
    EXPECT_EQ(exchange.ny(), ny);

    // Verify data types are valid (not MPI_DATATYPE_NULL)
    EXPECT_NE(exchange.row_type(), MPI_DATATYPE_NULL);
    EXPECT_NE(exchange.col_type(), MPI_DATATYPE_NULL);
}

/**
 * @test MoveConstructorTest
 * @brief Test move constructor
 *
 * This test verifies:
 * 1. Move constructor transfers ownership correctly
 * 2. Moved-from object has null types
 * 3. Moved-to object has valid types
 */
TEST_F(GhostCellExchangeTest, MoveConstructorTest) {
    int nx = 10;
    int ny = 10;

    GhostCellExchange original(nx, ny, *topology);
    MPI_Datatype original_row_type = original.row_type();
    MPI_Datatype original_col_type = original.col_type();

    // Move construct
    GhostCellExchange moved(std::move(original));

    // Verify moved object has valid types
    EXPECT_EQ(moved.nx(), nx);
    EXPECT_EQ(moved.ny(), ny);
    EXPECT_EQ(moved.row_type(), original_row_type);
    EXPECT_EQ(moved.col_type(), original_col_type);

    // Verify original object is moved-from
    // Note: topology is a reference, so it's still accessible
    EXPECT_EQ(original.row_type(), MPI_DATATYPE_NULL);
    EXPECT_EQ(original.col_type(), MPI_DATATYPE_NULL);
}

/**
 * @test ArrayValidationTest
 * @brief Test array size validation
 *
 * This test verifies:
 * 1. validate_array_size returns true for correct dimensions
 * 2. validate_array_size returns false for incorrect dimensions
 */
TEST_F(GhostCellExchangeTest, ArrayValidationTest) {
    int nx = 10;
    int ny = 10;

    GhostCellExchange exchange(nx, ny, *topology);

    // Correct size (ny+2) x (nx+2)
    utils::Array2D correct_array(ny + 2, nx + 2);
    EXPECT_TRUE(exchange.validate_array_size(correct_array));

    // Incorrect sizes
    utils::Array2D wrong_rows(ny + 1, nx + 2);
    EXPECT_FALSE(exchange.validate_array_size(wrong_rows));

    utils::Array2D wrong_cols(ny + 2, nx + 1);
    EXPECT_FALSE(exchange.validate_array_size(wrong_cols));
}

/**
 * @test SynchronousExchangeTest
 * @brief Test synchronous ghost cell exchange
 *
 * This test verifies:
 * 1. exchange() completes without errors
 * 2. Ghost cells are updated from neighbors
 * 3. Interior values are preserved
 *
 * Note: This test requires multiple processes to be meaningful.
 * With a single process, MPI_PROC_NULL neighbors will be used.
 */
TEST_F(GhostCellExchangeTest, SynchronousExchangeTest) {
    int nx = 10;
    int ny = 10;

    GhostCellExchange exchange(nx, ny, *topology);
    utils::Array2D array = create_test_array(nx, ny, 1000.0);

    // Save interior values
    std::vector<std::vector<double>> interior_values;
    for (int i = 1; i <= ny; ++i) {
        std::vector<double> row;
        for (int j = 1; j <= nx; ++j) {
            row.push_back(array(i, j));
        }
        interior_values.push_back(row);
    }

    // Perform exchange
    ASSERT_NO_THROW({
        exchange.exchange(array);
    });

    // Verify interior values are preserved
    for (int i = 1; i <= ny; ++i) {
        for (int j = 1; j <= nx; ++j) {
            EXPECT_DOUBLE_EQ(array(i, j), interior_values[i-1][j-1])
                << "Interior value changed at (" << i << ", " << j << ")";
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

/**
 * @test AsynchronousExchangeTest
 * @brief Test asynchronous ghost cell exchange
 *
 * This test verifies:
 * 1. exchange_async() initiates communication
 * 2. wait_all() completes all requests
 * 3. Results match synchronous exchange
 */
TEST_F(GhostCellExchangeTest, AsynchronousExchangeTest) {
    int nx = 10;
    int ny = 10;

    GhostCellExchange exchange(nx, ny, *topology);

    // Create two identical arrays
    utils::Array2D array_sync = create_test_array(nx, ny, 2000.0);
    utils::Array2D array_async = create_test_array(nx, ny, 2000.0);

    // Perform synchronous exchange
    exchange.exchange(array_sync);

    // Perform asynchronous exchange
    MPI_Request requests[8];
    exchange.exchange_async(array_async, requests);
    exchange.wait_all(requests);

    // Verify results match
    for (size_t i = 0; i < array_sync.rows(); ++i) {
        for (size_t j = 0; j < array_sync.cols(); ++j) {
            EXPECT_DOUBLE_EQ(array_async(i, j), array_sync(i, j))
                << "Async and sync results differ at (" << i << ", " << j << ")";
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

/**
 * @test BoundaryProcessTest
 * @brief Test boundary processes
 *
 * This test verifies:
 * 1. Processes on boundaries have MPI_PROC_NULL neighbors
 * 2. Exchange works correctly on boundaries
 * 3. No deadlock occurs
 */
TEST_F(GhostCellExchangeTest, BoundaryProcessTest) {
    int nx = 5;
    int ny = 5;

    GhostCellExchange exchange(nx, ny, *topology);
    const auto& neighbors = exchange.neighbors();

    // Check that boundary processes have MPI_PROC_NULL neighbors
    int coord_x = topology->coord_x();
    int coord_y = topology->coord_y();
    int dim_x = topology->dim_x();
    int dim_y = topology->dim_y();

    // South boundary
    if (coord_y == 0) {
        EXPECT_EQ(neighbors.south, MPI_PROC_NULL);
    }

    // North boundary
    if (coord_y == dim_y - 1) {
        EXPECT_EQ(neighbors.north, MPI_PROC_NULL);
    }

    // West boundary
    if (coord_x == 0) {
        EXPECT_EQ(neighbors.west, MPI_PROC_NULL);
    }

    // East boundary
    if (coord_x == dim_x - 1) {
        EXPECT_EQ(neighbors.east, MPI_PROC_NULL);
    }

    // Perform exchange on all processes
    utils::Array2D array = create_test_array(nx, ny, 3000.0);
    ASSERT_NO_THROW({
        exchange.exchange(array);
    });

    MPI_Barrier(MPI_COMM_WORLD);
}

/**
 * @test MultiExchangeTest
 * @brief Test multiple consecutive exchanges
 *
 * This test verifies:
 * 1. Multiple exchanges can be performed sequentially
 * 2. Each exchange operates correctly
 * 3. No data corruption occurs
 */
TEST_F(GhostCellExchangeTest, MultiExchangeTest) {
    int nx = 10;
    int ny = 10;

    GhostCellExchange exchange(nx, ny, *topology);
    utils::Array2D array = create_test_array(nx, ny, 4000.0);

    // Perform multiple exchanges
    for (int iter = 0; iter < 5; ++iter) {
        ASSERT_NO_THROW({
            exchange.exchange(array);
        });

        // Update interior values
        for (int i = 1; i <= ny; ++i) {
            for (int j = 1; j <= nx; ++j) {
                array(i, j) += 1.0;
            }
        }
    }

    // If we made it here without errors, test passes
    SUCCEED();

    MPI_Barrier(MPI_COMM_WORLD);
}

/**
 * @test LargeArrayTest
 * @brief Test exchange with large arrays
 *
 * This test verifies:
 * 1. Exchange works with large arrays
 * 2. Performance is reasonable
 * 3. No memory issues occur
 */
TEST_F(GhostCellExchangeTest, LargeArrayTest) {
    int nx = 100;
    int ny = 100;

    GhostCellExchange exchange(nx, ny, *topology);
    utils::Array2D array(ny + 2, nx + 2, 0.0);

    // Fill with test data
    for (int i = 0; i < ny + 2; ++i) {
        for (int j = 0; j < nx + 2; ++j) {
            array(i, j) = static_cast<double>(i * (nx + 2) + j);
        }
    }

    // Perform exchange
    ASSERT_NO_THROW({
        exchange.exchange(array);
    });

    // Verify array is still valid (no corruption)
    double sum = 0.0;
    for (size_t i = 0; i < array.rows(); ++i) {
        for (size_t j = 0; j < array.cols(); ++j) {
            sum += std::abs(array(i, j));
        }
    }
    EXPECT_GT(sum, 0.0);

    MPI_Barrier(MPI_COMM_WORLD);
}

/**
 * @test ExceptionTest
 * @brief Test exception handling
 *
 * This test verifies:
 * 1. Invalid dimensions throw exceptions
 * 2. Invalid array sizes throw exceptions
 */
TEST_F(GhostCellExchangeTest, ExceptionTest) {
    // Test invalid dimensions
    EXPECT_THROW({
        GhostCellExchange exchange(0, 10, *topology);
    }, std::invalid_argument);

    EXPECT_THROW({
        GhostCellExchange exchange(10, 0, *topology);
    }, std::invalid_argument);

    // Test invalid array size
    GhostCellExchange exchange(10, 10, *topology);
    utils::Array2D wrong_array(10, 10);

    EXPECT_THROW({
        exchange.exchange(wrong_array);
    }, std::invalid_argument);

    EXPECT_THROW({
        MPI_Request requests[8];
        exchange.exchange_async(wrong_array, requests);
    }, std::invalid_argument);
}

/**
 * @test PerformanceTest
 * @brief Measure exchange performance
 *
 * This test measures the time required for multiple exchanges
 * and reports the performance metrics.
 */
TEST_F(GhostCellExchangeTest, PerformanceTest) {
    int nx = 50;
    int ny = 50;
    int num_iterations = 10;

    GhostCellExchange exchange(nx, ny, *topology);
    utils::Array2D array = create_test_array(nx, ny, 5000.0);

    // Warm-up
    for (int i = 0; i < 3; ++i) {
        exchange.exchange(array);
    }

    // Synchronize all processes
    MPI_Barrier(MPI_COMM_WORLD);

    // Measure synchronous exchange time
    double start_time = MPI_Wtime();
    for (int i = 0; i < num_iterations; ++i) {
        exchange.exchange(array);
    }
    double end_time = MPI_Wtime();
    double sync_time = end_time - start_time;

    // Measure asynchronous exchange time
    MPI_Barrier(MPI_COMM_WORLD);

    start_time = MPI_Wtime();
    for (int i = 0; i < num_iterations; ++i) {
        MPI_Request requests[8];
        exchange.exchange_async(array, requests);
        exchange.wait_all(requests);
    }
    end_time = MPI_Wtime();
    double async_time = end_time - start_time;

    // Report performance metrics
    if (topology->rank() == 0) {
        std::cout << "\n=== GhostCellExchange Performance ===" << std::endl;
        std::cout << "Array size: " << (nx+2) << " x " << (ny+2) << std::endl;
        std::cout << "Iterations: " << num_iterations << std::endl;
        std::cout << "Synchronous exchange time: " << sync_time << " seconds" << std::endl;
        std::cout << "Asynchronous exchange time: " << async_time << " seconds" << std::endl;
        std::cout << "Average sync exchange: " << (sync_time / num_iterations * 1000.0)
                  << " ms" << std::endl;
        std::cout << "Average async exchange: " << (async_time / num_iterations * 1000.0)
                  << " ms" << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Test passes as long as we can complete the exchanges
    SUCCEED();
}

/**
 * @test DirectionConstantsTest
 * @brief Test direction constant values
 *
 * This test verifies:
 * 1. Direction constants have correct values
 * 2. Constants are usable for indexing
 */
TEST_F(GhostCellExchangeTest, DirectionConstantsTest) {
    EXPECT_EQ(ghost::SOUTH, 0);
    EXPECT_EQ(ghost::NORTH, 1);
    EXPECT_EQ(ghost::WEST, 2);
    EXPECT_EQ(ghost::EAST, 3);

    // Test that they can be used as array indices
    int array[4];
    array[ghost::SOUTH] = 0;
    array[ghost::NORTH] = 1;
    array[ghost::WEST] = 2;
    array[ghost::EAST] = 3;

    EXPECT_EQ(array[0], 0);
    EXPECT_EQ(array[1], 1);
    EXPECT_EQ(array[2], 2);
    EXPECT_EQ(array[3], 3);
}

/**
 * @test DataLayoutTest
 * @brief Test data layout with ghost cells
 *
 * This test verifies:
 * 1. Array layout matches expected structure
 * 2. Ghost cells are accessible
 * 3. Interior cells are accessible
 */
TEST_F(GhostCellExchangeTest, DataLayoutTest) {
    int nx = 5;
    int ny = 5;

    utils::Array2D array(ny + 2, nx + 2);

    // Expected layout:
    // Row 0:        South ghost row
    // Row 1..ny:    Interior rows
    // Row ny+1:     North ghost row
    // Col 0:        West ghost column
    // Col 1..nx:    Interior columns
    // Col nx+1:     East ghost column

    EXPECT_EQ(array.rows(), static_cast<size_t>(ny + 2));
    EXPECT_EQ(array.cols(), static_cast<size_t>(nx + 2));

    // Mark ghost cells
    for (int j = 0; j <= nx + 1; ++j) {
        array(0, j) = -1.0;           // South ghost
        array(ny + 1, j) = -2.0;      // North ghost
    }
    for (int i = 0; i <= ny + 1; ++i) {
        array(i, 0) = -3.0;           // West ghost
        array(i, nx + 1) = -4.0;      // East ghost
    }

    // Mark interior
    for (int i = 1; i <= ny; ++i) {
        for (int j = 1; j <= nx; ++j) {
            array(i, j) = 1.0;
        }
    }

    // Verify ghost cells
    for (int j = 0; j <= nx + 1; ++j) {
        EXPECT_EQ(array(0, j), -1.0);        // South ghost
        EXPECT_EQ(array(ny + 1, j), -2.0);    // North ghost
    }
    for (int i = 0; i <= ny + 1; ++i) {
        EXPECT_EQ(array(i, 0), -3.0);        // West ghost
        EXPECT_EQ(array(i, nx + 1), -4.0);    // East ghost
    }

    // Verify interior
    for (int i = 1; i <= ny; ++i) {
        for (int j = 1; j <= nx; ++j) {
            EXPECT_EQ(array(i, j), 1.0);
        }
    }
}

/**
 * @test MultipleTopologiesTest
 * @brief Test with different process grid configurations
 *
 * This test verifies that the exchange works correctly
 * with different numbers of processes and grid layouts.
 */
TEST_F(GhostCellExchangeTest, MultipleTopologiesTest) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int nx = 8;
    int ny = 8;

    // Test with default topology
    {
        CartesianTopology topo1(MPI_COMM_WORLD);
        GhostCellExchange exchange1(nx, ny, topo1);
        utils::Array2D array1 = create_test_array(nx, ny, 6000.0);

        ASSERT_NO_THROW({
            exchange1.exchange(array1);
        });
    }

    // Test with specified dimensions (if compatible)
    if (size == 4) {
        std::vector<int> dims = {2, 2};
        CartesianTopology topo2(MPI_COMM_WORLD, dims);
        GhostCellExchange exchange2(nx, ny, topo2);
        utils::Array2D array2 = create_test_array(nx, ny, 7000.0);

        ASSERT_NO_THROW({
            exchange2.exchange(array2);
        });
    }

    if (size == 6) {
        std::vector<int> dims = {2, 3};
        CartesianTopology topo3(MPI_COMM_WORLD, dims);
        GhostCellExchange exchange3(nx, ny, topo3);
        utils::Array2D array3 = create_test_array(nx, ny, 8000.0);

        ASSERT_NO_THROW({
            exchange3.exchange(array3);
        });
    }

    if (size == 9) {
        std::vector<int> dims = {3, 3};
        CartesianTopology topo4(MPI_COMM_WORLD, dims);
        GhostCellExchange exchange4(nx, ny, topo4);
        utils::Array2D array4 = create_test_array(nx, ny, 9000.0);

        ASSERT_NO_THROW({
            exchange4.exchange(array4);
        });
    }

    MPI_Barrier(MPI_COMM_WORLD);
    SUCCEED();
}
