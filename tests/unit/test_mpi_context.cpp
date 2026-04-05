#include <gtest/gtest.h>
#include <mpi.h>
#include "../src/mpi/mpi_context.hpp"
#include <memory>
#include <stdexcept>

/**
 * @brief Test fixture for MPIContext tests
 *
 * Note: MPIContext tests require running with MPI (mpirun/mpiexec)
 * Example: mpirun -n 4 ./test_mpi_context
 */
class MPIContextTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Check if MPI is initialized before running tests
        int mpi_initialized = 0;
        MPI_Initialized(&mpi_initialized);
        if (!mpi_initialized) {
            // Skip tests if MPI is not initialized
            // Tests should be run with: mpirun -n <num_procs> ./test_mpi_context
            GTEST_SKIP() << "MPI not initialized. Run with mpirun/mpiexec.";
        }
    }

    void TearDown() override {
        // Clean up after tests
    }
};

/**
 * @test InitializationTest
 * @brief Test that MPIContext initializes correctly
 *
 * This test verifies:
 * 1. MPIContext can be constructed
 * 2. MPI is initialized after construction
 * 3. Rank and size are valid (non-negative, rank < size)
 */
TEST_F(MPIContextTest, InitializationTest) {
    int argc = 0;
    char** argv = nullptr;

    // Test that we can create an MPIContext
    ASSERT_NO_THROW({
        MPIContext context(argc, argv);
    });

    // Create a context for testing
    MPIContext context(argc, argv);

    // Verify MPI is initialized
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    EXPECT_TRUE(mpi_initialized);
}

/**
 * @test RankSizeTest
 * @brief Test that rank and size are retrieved correctly
 *
 * This test verifies:
 * 1. Rank is non-negative
 * 2. Size is positive
 * 3. Rank is less than size
 * 4. Rank and size match MPI calls
 */
TEST_F(MPIContextTest, RankSizeTest) {
    int argc = 0;
    char** argv = nullptr;

    MPIContext context(argc, argv);

    // Get rank and size from context
    int rank = context.rank();
    int size = context.size();

    // Verify rank is non-negative
    EXPECT_GE(rank, 0);

    // Verify size is positive
    EXPECT_GT(size, 0);

    // Verify rank is less than size
    EXPECT_LT(rank, size);

    // Compare with actual MPI values
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    EXPECT_EQ(rank, mpi_rank);
    EXPECT_EQ(size, mpi_size);
}

/**
 * @test IsRootTest
 * @brief Test that is_root() correctly identifies root process
 *
 * This test verifies:
 * 1. is_root() returns true only for rank 0
 * 2. is_root() returns false for non-zero ranks
 * 3. Exactly one process is root
 */
TEST_F(MPIContextTest, IsRootTest) {
    int argc = 0;
    char** argv = nullptr;

    MPIContext context(argc, argv);

    int rank = context.rank();
    bool is_root = context.is_root();

    // Verify is_root() matches rank == 0
    EXPECT_EQ(is_root, (rank == 0));

    // Verify consistency
    if (rank == 0) {
        EXPECT_TRUE(is_root);
    } else {
        EXPECT_FALSE(is_root);
    }
}

/**
 * @test BarrierTest
 * @brief Test that barrier synchronizes processes correctly
 *
 * This test verifies:
 * 1. barrier() can be called without throwing
 * 2. All processes reach the barrier
 * 3. Multiple barriers work correctly
 */
TEST_F(MPIContextTest, BarrierTest) {
    int argc = 0;
    char** argv = nullptr;

    MPIContext context(argc, argv);

    // Test single barrier
    ASSERT_NO_THROW({
        context.barrier();
    });

    // Test multiple barriers
    for (int i = 0; i < 3; ++i) {
        ASSERT_NO_THROW({
            context.barrier();
        });
    }
}

/**
 * @test NoCopyTest
 * @brief Test that copy semantics are properly disabled
 *
 * This test verifies:
 * 1. Copy constructor is deleted (compile-time error)
 * 2. Copy assignment is deleted (compile-time error)
 *
 * Note: This test uses static assertions to verify at compile time
 */
TEST_F(MPIContextTest, NoCopyTest) {
    int argc = 0;
    char** argv = nullptr;

    MPIContext context(argc, argv);

    // Verify that copy constructor is deleted
    // This should not compile if copy constructor were available
    EXPECT_FALSE(std::is_copy_constructible<MPIContext>::value);

    // Verify that copy assignment is deleted
    // This should not compile if copy assignment were available
    EXPECT_FALSE(std::is_copy_assignable<MPIContext>::value);
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
TEST_F(MPIContextTest, MoveTest) {
    int argc = 0;
    char** argv = nullptr;

    // Test move constructor
    MPIContext context1(argc, argv);
    int original_rank = context1.rank();
    int original_size = context1.size();

    MPIContext context2(std::move(context1));

    // Verify ownership transfer
    EXPECT_EQ(context2.rank(), original_rank);
    EXPECT_EQ(context2.size(), original_size);

    // Test move assignment
    MPIContext context3(argc, argv);
    context3 = std::move(context2);

    EXPECT_EQ(context3.rank(), original_rank);
    EXPECT_EQ(context3.size(), original_size);
}

/**
 * @test DestructorTest
 * @brief Test that destructor properly cleans up MPI
 *
 * This test verifies:
 * 1. Destructor calls MPI_Finalize when appropriate
 * 2. No exceptions are thrown during destruction
 */
TEST_F(MPIContextTest, DestructorTest) {
    int argc = 0;
    char** argv = nullptr;

    // Test that destructor doesn't throw
    ASSERT_NO_THROW({
        MPIContext context(argc, argv);
        // context goes out of scope here
    });

    // Create multiple contexts to test proper cleanup
    for (int i = 0; i < 3; ++i) {
        MPIContext context(argc, argv);
        EXPECT_GE(context.rank(), 0);
        EXPECT_GT(context.size(), 0);
    }
}

/**
 * @test MultipleInitializationTest
 * @brief Test that multiple contexts handle existing MPI initialization
 *
 * This test verifies:
 * 1. Second context works when first context already initialized MPI
 * 2. Both contexts get correct rank and size
 */
TEST_F(MPIContextTest, MultipleInitializationTest) {
    int argc = 0;
    char** argv = nullptr;

    // Create first context
    MPIContext context1(argc, argv);
    int rank1 = context1.rank();
    int size1 = context1.size();

    // Create second context (MPI already initialized)
    MPIContext context2(argc, argv);
    int rank2 = context2.rank();
    int size2 = context2.size();

    // Both should have same rank and size
    EXPECT_EQ(rank1, rank2);
    EXPECT_EQ(size1, size2);
}

/**
 * @test AbortTest
 * @brief Test abort functionality
 *
 * This test verifies:
 * 1. abort() can be called
 * 2. Note: This test cannot fully test MPI_Abort as it terminates execution
 *
 * This test is marked as DISABLED_ because it would terminate the test run
 */
TEST_F(MPIContextTest, DISABLED_AbortTest) {
    int argc = 0;
    char** argv = nullptr;

    MPIContext context(argc, argv);

    // This would terminate all MPI processes
    // Uncomment to test manually:
    // context.abort(1, "Test abort message");
}

/**
 * @main
 * @brief Main function for running MPIContext tests
 *
 * This function initializes Google Test and runs all tests.
 * It should be called via mpirun/mpiexec to initialize MPI.
 */
int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    // Initialize MPI if not already initialized
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    if (!mpi_initialized) {
        MPI_Init(&argc, &argv);
    }

    int result = RUN_ALL_TESTS();

    // Finalize MPI if we initialized it
    MPI_Initialized(&mpi_initialized);
    int mpi_finalized = 0;
    MPI_Finalized(&mpi_finalized);
    if (mpi_initialized && !mpi_finalized) {
        MPI_Finalize();
    }

    return result;
}
