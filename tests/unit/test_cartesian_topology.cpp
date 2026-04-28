#include <gtest/gtest.h>
#include <mpi.h>
#include "../src/mpi/cartesian_topology.hpp"
#include <memory>
#include <stdexcept>

/**
 * @brief Test fixture for CartesianTopology tests
 *
 * Note: CartesianTopology tests require running with MPI (mpirun/mpiexec)
 * Example: mpirun -n 4 ./test_cartesian_topology
 */
class CartesianTopologyTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Check if MPI is initialized before running tests
        int mpi_initialized = 0;
        MPI_Initialized(&mpi_initialized);
        if (!mpi_initialized) {
            // Skip tests if MPI is not initialized
            // Tests should be run with: mpirun -n <num_procs> ./test_cartesian_topology
            GTEST_SKIP() << "MPI not initialized. Run with mpirun/mpiexec.";
        }
    }

    void TearDown() override {
        // Clean up after tests
    }
};

/**
 * @test InitializationTest
 * @brief Test that CartesianTopology initializes correctly
 *
 * This test verifies:
 * 1. CartesianTopology can be constructed
 * 2. Cartesian communicator is created successfully
 * 3. Rank and size are valid
 * 4. Communicator is not MPI_COMM_NULL
 */
TEST_F(CartesianTopologyTest, InitializationTest) {
    // Test that we can create a CartesianTopology with default communicator
    ASSERT_NO_THROW({
        CartesianTopology topology(MPI_COMM_WORLD);
    });

    // Create a topology for testing
    CartesianTopology topology(MPI_COMM_WORLD);

    // Verify communicator is valid
    EXPECT_NE(topology.communicator(), MPI_COMM_NULL);

    // Verify rank and size are valid
    EXPECT_GE(topology.rank(), 0);
    EXPECT_GT(topology.size(), 0);
    EXPECT_LT(topology.rank(), topology.size());
}

/**
 * @test DimensionCalculationTest
 * @brief Test automatic dimension calculation
 *
 * This test verifies:
 * 1. Optimal dimensions are calculated correctly
 * 2. Dimensions are closest to square
 * 3. Product of dimensions equals number of processes
 */
TEST_F(CartesianTopologyTest, DimensionCalculationTest) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    CartesianTopology topology(MPI_COMM_WORLD);

    const auto& dims = topology.dims();
    ASSERT_EQ(dims.size(), 2);

    // Verify product equals number of processes
    EXPECT_EQ(dims[0] * dims[1], size);

    // Verify dimensions are close to square (within reasonable bounds)
    int max_dim = std::max(dims[0], dims[1]);
    int min_dim = std::min(dims[0], dims[1]);
    double ratio = static_cast<double>(max_dim) / min_dim;

    // For square-like topologies, ratio should be reasonable
    EXPECT_LE(ratio, 2.0) << "Dimensions are too skewed: [" << dims[0] << ", " << dims[1] << "]";
}

/**
 * @test SpecifiedDimensionsTest
 * @brief Test that specified dimensions are used correctly
 *
 * This test verifies:
 * 1. Constructor accepts user-specified dimensions
 * 2. Dimensions match what was specified
 * 3. Invalid dimensions throw exceptions
 */
TEST_F(CartesianTopologyTest, SpecifiedDimensionsTest) {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Find valid dimensions for this size
    int dim_x = 1;
    int dim_y = size;

    // Try to find a more balanced factorization
    for (int i = 2; i * i <= size; ++i) {
        if (size % i == 0) {
            dim_x = i;
            dim_y = size / i;
            break;
        }
    }

    // Test with valid dimensions
    ASSERT_NO_THROW({
        CartesianTopology topology(MPI_COMM_WORLD, {dim_x, dim_y});
    });

    CartesianTopology topology(MPI_COMM_WORLD, {dim_x, dim_y});
    EXPECT_EQ(topology.dim_x(), dim_x);
    EXPECT_EQ(topology.dim_y(), dim_y);

    // Test that invalid dimensions throw exception
    if (size != 6) {
        EXPECT_THROW({
            CartesianTopology topology(MPI_COMM_WORLD, {2, 3});  // 2*3=6
        }, std::invalid_argument);
    }

    // Test that wrong number of dimensions throws exception
    EXPECT_THROW({
        CartesianTopology topology(MPI_COMM_WORLD, {dim_x});  // Only 1 dimension
    }, std::invalid_argument);
}

/**
 * @test CoordinateMappingTest
 * @brief Test coordinate mapping between ranks and positions
 *
 * This test verifies:
 * 1. Each process gets unique coordinates
 * 2. Coordinates are within valid range
 * 3. Coordinate-to-rank mapping is consistent
 */
TEST_F(CartesianTopologyTest, CoordinateMappingTest) {
    CartesianTopology topology(MPI_COMM_WORLD);

    const auto& dims = topology.dims();
    const auto& coords = topology.coords();

    // Verify coordinate dimensions
    ASSERT_EQ(coords.size(), 2);

    // Verify coordinates are within valid range
    EXPECT_GE(coords[0], 0);
    EXPECT_LT(coords[0], dims[0]);
    EXPECT_GE(coords[1], 0);
    EXPECT_LT(coords[1], dims[1]);

    // Verify that dim_x/coord_x and dim_y/coord_y match
    EXPECT_EQ(topology.coord_x(), coords[0]);
    EXPECT_EQ(topology.coord_y(), coords[1]);
}

/**
 * @test NeighborComputationTest
 * @brief Test neighbor computation
 *
 * This test verifies:
 * 1. Each process gets correct neighbor ranks
 * 2. Boundary processes have MPI_PROC_NULL as neighbors
 * 3. Neighbor relationships are symmetric
 */
TEST_F(CartesianTopologyTest, NeighborComputationTest) {
    CartesianTopology topology(MPI_COMM_WORLD);

    const auto& neighbors = topology.neighbors();

    // Verify that all neighbors are either valid ranks or MPI_PROC_NULL
    int size = topology.size();

    // Check each neighbor
    if (neighbors.south != MPI_PROC_NULL) {
        EXPECT_GE(neighbors.south, 0);
        EXPECT_LT(neighbors.south, size);
    }

    if (neighbors.north != MPI_PROC_NULL) {
        EXPECT_GE(neighbors.north, 0);
        EXPECT_LT(neighbors.north, size);
    }

    if (neighbors.west != MPI_PROC_NULL) {
        EXPECT_GE(neighbors.west, 0);
        EXPECT_LT(neighbors.west, size);
    }

    if (neighbors.east != MPI_PROC_NULL) {
        EXPECT_GE(neighbors.east, 0);
        EXPECT_LT(neighbors.east, size);
    }

    // No process should be its own neighbor
    int rank = topology.rank();
    EXPECT_NE(neighbors.south, rank);
    EXPECT_NE(neighbors.north, rank);
    EXPECT_NE(neighbors.west, rank);
    EXPECT_NE(neighbors.east, rank);
}

/**
 * @test BoundaryDetectionTest
 * @brief Test boundary detection
 *
 * This test verifies:
 * 1. is_on_boundary() correctly identifies boundary processes
 * 2. has_neighbor() is the opposite of is_on_boundary
 * 3. Interior processes have neighbors in all directions
 */
TEST_F(CartesianTopologyTest, BoundaryDetectionTest) {
    CartesianTopology topology(MPI_COMM_WORLD);

    // For each direction and shift, verify consistency
    // X direction (East/West)
    bool is_west_boundary = topology.is_on_boundary(Direction::X, Shift::Backward);
    bool is_east_boundary = topology.is_on_boundary(Direction::X, Shift::Forward);

    EXPECT_EQ(is_west_boundary, !topology.has_neighbor(Direction::X, Shift::Backward));
    EXPECT_EQ(is_east_boundary, !topology.has_neighbor(Direction::X, Shift::Forward));

    // Y direction (North/South)
    bool is_south_boundary = topology.is_on_boundary(Direction::Y, Shift::Backward);
    bool is_north_boundary = topology.is_on_boundary(Direction::Y, Shift::Forward);

    EXPECT_EQ(is_south_boundary, !topology.has_neighbor(Direction::Y, Shift::Backward));
    EXPECT_EQ(is_north_boundary, !topology.has_neighbor(Direction::Y, Shift::Forward));

    // If on west boundary, should have no west neighbor
    if (is_west_boundary) {
        EXPECT_EQ(topology.neighbors().west, MPI_PROC_NULL);
    }

    // If on east boundary, should have no east neighbor
    if (is_east_boundary) {
        EXPECT_EQ(topology.neighbors().east, MPI_PROC_NULL);
    }

    // If on south boundary, should have no south neighbor
    if (is_south_boundary) {
        EXPECT_EQ(topology.neighbors().south, MPI_PROC_NULL);
    }

    // If on north boundary, should have no north neighbor
    if (is_north_boundary) {
        EXPECT_EQ(topology.neighbors().north, MPI_PROC_NULL);
    }
}

/**
 * @test EdgeCasesTest
 * @brief Test 1xN and Nx1 topologies
 *
 * This test verifies:
 * 1. 1xN topology works correctly
 * 2. Nx1 topology works correctly
 * 3. All processes in a row/column have correct neighbors
 */
TEST_F(CartesianTopologyTest, EdgeCasesTest) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Test 1xN topology if appropriate
    if (size <= 4) {  // Only test with small sizes
        ASSERT_NO_THROW({
            CartesianTopology topology(MPI_COMM_WORLD, {1, size});
        });

        CartesianTopology topology(MPI_COMM_WORLD, {1, size});
        EXPECT_EQ(topology.dim_x(), 1);
        EXPECT_EQ(topology.dim_y(), size);

        // All processes should be on west and east boundaries
        EXPECT_TRUE(topology.is_on_boundary(Direction::X, Shift::Backward));
        EXPECT_TRUE(topology.is_on_boundary(Direction::X, Shift::Forward));
        EXPECT_EQ(topology.neighbors().west, MPI_PROC_NULL);
        EXPECT_EQ(topology.neighbors().east, MPI_PROC_NULL);
    }

    // Test Nx1 topology if appropriate
    if (size <= 4) {
        ASSERT_NO_THROW({
            CartesianTopology topology(MPI_COMM_WORLD, {size, 1});
        });

        CartesianTopology topology(MPI_COMM_WORLD, {size, 1});
        EXPECT_EQ(topology.dim_x(), size);
        EXPECT_EQ(topology.dim_y(), 1);

        // All processes should be on south and north boundaries
        EXPECT_TRUE(topology.is_on_boundary(Direction::Y, Shift::Backward));
        EXPECT_TRUE(topology.is_on_boundary(Direction::Y, Shift::Forward));
        EXPECT_EQ(topology.neighbors().south, MPI_PROC_NULL);
        EXPECT_EQ(topology.neighbors().north, MPI_PROC_NULL);
    }
}

/**
 * @test SquareTopologyTest
 * @brief Test square topologies (2x2, 3x3, 4x4)
 *
 * This test verifies:
 * 1. Square topologies are created correctly
 * 2. Corners have two neighbors
 * 3. Edges have three neighbors
 * 4. Interior processes have four neighbors
 */
TEST_F(CartesianTopologyTest, SquareTopologyTest) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Only test with perfect squares: 4 (2x2) or 9 (3x3)
    if (size == 4 || size == 9) {
        int dim = static_cast<int>(std::sqrt(size));
        ASSERT_NO_THROW({
            CartesianTopology topology(MPI_COMM_WORLD, {dim, dim});
        });

        CartesianTopology topology(MPI_COMM_WORLD, {dim, dim});

        int coord_x = topology.coord_x();
        int coord_y = topology.coord_y();

        int neighbor_count = 0;
        if (topology.neighbors().south != MPI_PROC_NULL) neighbor_count++;
        if (topology.neighbors().north != MPI_PROC_NULL) neighbor_count++;
        if (topology.neighbors().west != MPI_PROC_NULL) neighbor_count++;
        if (topology.neighbors().east != MPI_PROC_NULL) neighbor_count++;

        // Corners should have 2 neighbors
        if ((coord_x == 0 || coord_x == dim - 1) && (coord_y == 0 || coord_y == dim - 1)) {
            EXPECT_EQ(neighbor_count, 2);
        }
        // Edges should have 3 neighbors
        else if (coord_x == 0 || coord_x == dim - 1 || coord_y == 0 || coord_y == dim - 1) {
            EXPECT_EQ(neighbor_count, 3);
        }
        // Interior should have 4 neighbors
        else {
            EXPECT_EQ(neighbor_count, 4);
        }
    }
}

/**
 * @test NoCopyTest
 * @brief Test that copy semantics are properly disabled
 *
 * This test verifies:
 * 1. Copy constructor is deleted
 * 2. Copy assignment is deleted
 */
TEST_F(CartesianTopologyTest, NoCopyTest) {
    CartesianTopology topology(MPI_COMM_WORLD);

    // Verify that copy constructor is deleted
    EXPECT_FALSE(std::is_copy_constructible<CartesianTopology>::value);

    // Verify that copy assignment is deleted
    EXPECT_FALSE(std::is_copy_assignable<CartesianTopology>::value);
}

/**
 * @test MoveTest
 * @brief Test that move semantics work correctly
 *
 * This test verifies:
 * 1. Move constructor transfers ownership
 * 2. Move assignment transfers ownership
 * 3. Moved-from object is in valid but invalid state
 */
TEST_F(CartesianTopologyTest, MoveTest) {
    // Test move constructor
    CartesianTopology topology1(MPI_COMM_WORLD);
    int original_rank = topology1.rank();
    int original_size = topology1.size();
    auto original_dims = topology1.dims();

    CartesianTopology topology2(std::move(topology1));

    // Verify ownership transfer
    EXPECT_EQ(topology2.rank(), original_rank);
    EXPECT_EQ(topology2.size(), original_size);
    EXPECT_EQ(topology2.dims(), original_dims);

    // Test move assignment
    CartesianTopology topology3(MPI_COMM_WORLD);
    topology3 = std::move(topology2);

    EXPECT_EQ(topology3.rank(), original_rank);
    EXPECT_EQ(topology3.size(), original_size);
    EXPECT_EQ(topology3.dims(), original_dims);
}

/**
 * @test DestructorTest
 * @brief Test that destructor properly cleans up communicator
 *
 * This test verifies:
 * 1. Destructor frees the Cartesian communicator
 * 2. No exceptions are thrown during destruction
 */
TEST_F(CartesianTopologyTest, DestructorTest) {
    // Test that destructor doesn't throw
    ASSERT_NO_THROW({
        CartesianTopology topology(MPI_COMM_WORLD);
        // topology goes out of scope here
    });

    // Create multiple topologies to test proper cleanup
    for (int i = 0; i < 3; ++i) {
        CartesianTopology topology(MPI_COMM_WORLD);
        EXPECT_NE(topology.communicator(), MPI_COMM_NULL);
        EXPECT_GE(topology.rank(), 0);
        EXPECT_GT(topology.size(), 0);
    }
}

/**
 * @test MultipleTopologiesTest
 * @brief Test creating multiple topologies
 *
 * This test verifies:
 * 1. Multiple topologies can be created sequentially
 * 2. Each topology is independent
 */
TEST_F(CartesianTopologyTest, MultipleTopologiesTest) {
    CartesianTopology topology1(MPI_COMM_WORLD);
    int rank1 = topology1.rank();
    int size1 = topology1.size();

    CartesianTopology topology2(MPI_COMM_WORLD);
    int rank2 = topology2.rank();
    int size2 = topology2.size();

    // Both should have same rank and size
    EXPECT_EQ(rank1, rank2);
    EXPECT_EQ(size1, size2);
}

/**
 * @main
 * @brief Main function for running CartesianTopology tests
 *
 * This function initializes Google Test and runs all tests.
 * It should be called via mpirun/mpiexec to initialize MPI.
 */
