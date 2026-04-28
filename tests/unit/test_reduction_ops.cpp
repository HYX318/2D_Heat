/**
 * @file test_reduction_ops.cpp
 * @brief Unit tests for ReductionOps class
 */

#include <gtest/gtest.h>
#include <mpi_context.hpp>
#include <reduction_ops.hpp>
#include <vector>
#include <cmath>
#include <numeric>
#include "../test_common.hpp"

using namespace test_utils;

// Test fixture for ReductionOps tests
class ReductionOpsTest : public ::testing::Test {
protected:
    void SetUp() override {
        skip_if_no_mpi();
        ops = std::make_unique<ReductionOps>();
    }

    void TearDown() override {
        ops.reset();
    }

    std::unique_ptr<ReductionOps> ops;
};

// ============================================================================
// Scalar Reduce Tests
// ============================================================================

TEST_F(ReductionOpsTest, ReduceSumInt) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    int result = ops->reduce(rank * 2, ReductionOp::SUM, 0);

    if (rank == 0) {
        int expected = 0 + 2 + 4 + 6;  // 0*2 + 1*2 + 2*2 + 3*2
        EXPECT_EQ(result, expected);
    }
}

TEST_F(ReductionOpsTest, ReduceMaxInt) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    int result = ops->reduce(rank, ReductionOp::MAX, 0);

    if (rank == 0) {
        EXPECT_EQ(result, 3);  // Maximum rank
    }
}

TEST_F(ReductionOpsTest, ReduceMinInt) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    int result = ops->reduce(rank, ReductionOp::MIN, 1);

    if (rank == 1) {
        EXPECT_EQ(result, 0);  // Minimum rank
    }
}

TEST_F(ReductionOpsTest, ReduceProdInt) {
    ASSERT_MPI_PROCS(3);

    int rank = get_mpi_rank();
    int value = rank + 1;  // 1, 2, 3
    int result = ops->reduce(value, ReductionOp::PROD, 0);

    if (rank == 0) {
        EXPECT_EQ(result, 6);  // 1 * 2 * 3
    }
}

TEST_F(ReductionOpsTest, ReduceSumDouble) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    double value = rank * 0.5;
    double result = ops->reduce(value, ReductionOp::SUM, 0);

    if (rank == 0) {
        double expected = 0.0 + 0.5 + 1.0 + 1.5;
        EXPECT_DOUBLE_EQ(result, expected);
    }
}

// ============================================================================
// Allreduce Tests
// ============================================================================

TEST_F(ReductionOpsTest, AllreduceSumInt) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    int result = ops->allreduce(rank, ReductionOp::SUM);

    int expected = 0 + 1 + 2 + 3;
    EXPECT_EQ(result, expected);
}

TEST_F(ReductionOpsTest, AllreduceMaxDouble) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    double value = rank * 1.5;
    double result = ops->allreduce(value, ReductionOp::MAX);

    EXPECT_DOUBLE_EQ(result, 4.5);  // 3 * 1.5
}

TEST_F(ReductionOpsTest, AllreduceMinFloat) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    float value = static_cast<float>(rank + 10);
    float result = ops->allreduce(value, ReductionOp::MIN);

    EXPECT_FLOAT_EQ(result, 10.0f);
}

TEST_F(ReductionOpsTest, AllreduceLogicalAnd) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    int value = (rank % 2 == 0) ? 1 : 0;  // 1, 0, 1, 0
    int result = ops->allreduce(value, ReductionOp::LAND);

    EXPECT_EQ(result, 0);  // Not all are true
}

TEST_F(ReductionOpsTest, AllreduceLogicalOr) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    int value = (rank == 1) ? 1 : 0;  // 0, 1, 0, 0
    int result = ops->allreduce(value, ReductionOp::LOR);

    EXPECT_EQ(result, 1);  // At least one is true
}

TEST_F(ReductionOpsTest, AllreduceBitwiseAnd) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    unsigned int value = 0xF0 | (rank & 0x0F);
    unsigned int result = ops->allreduce(value, ReductionOp::BAND);

    EXPECT_EQ(result, 0xF0);  // Common bits
}

TEST_F(ReductionOpsTest, AllreduceBitwiseOr) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    unsigned int value = rank;
    unsigned int result = ops->allreduce(value, ReductionOp::BOR);

    EXPECT_EQ(result, 0x03);  // 0 | 1 | 2 | 3 = 3
}

// ============================================================================
// Scan Tests
// ============================================================================

TEST_F(ReductionOpsTest, ScanSumInclusive) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    int result = ops->scan(rank + 1, ReductionOp::SUM);

    // Expected: 1, 3, 6, 10 for ranks 0, 1, 2, 3
    int expected = (rank + 1) * (rank + 2) / 2;
    EXPECT_EQ(result, expected);
}

TEST_F(ReductionOpsTest, ScanMaxInclusive) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    int result = ops->scan(rank, ReductionOp::MAX);

    // Expected: 0, 1, 2, 3 for ranks 0, 1, 2, 3
    EXPECT_EQ(result, rank);
}

TEST_F(ReductionOpsTest, ExscanSumExclusive) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    int result = ops->exscan(rank + 1, ReductionOp::SUM);

    // Expected: 0, 1, 3, 6 for ranks 0, 1, 2, 3
    int expected = rank * (rank + 1) / 2;
    EXPECT_EQ(result, expected);
}

TEST_F(ReductionOpsTest, ExscanMaxExclusive) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    int result = ops->exscan(rank, ReductionOp::MAX);

    // Expected: min value for rank 0, 0, 1, 2 for ranks 0, 1, 2, 3
    int expected = (rank == 0) ? std::numeric_limits<int>::lowest() : rank - 1;
    EXPECT_EQ(result, expected);
}

// ============================================================================
// Array Reduce Tests
// ============================================================================

TEST_F(ReductionOpsTest, ReduceArraySum) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    std::vector<int> values = {rank, rank * 2, rank * 3};

    std::vector<int> result = ops->reduce_array(values, ReductionOp::SUM, 0);

    if (rank == 0) {
        std::vector<int> expected = {6, 12, 18};  // Sum of each column
        EXPECT_EQ(result.size(), expected.size());
        for (size_t i = 0; i < result.size(); ++i) {
            EXPECT_EQ(result[i], expected[i]);
        }
    } else {
        EXPECT_TRUE(result.empty());
    }
}

TEST_F(ReductionOpsTest, AllreduceArrayMax) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    std::vector<double> values = {rank * 1.0, rank * 2.0, rank * 3.0};

    std::vector<double> result = ops->allreduce_array(values, ReductionOp::MAX);

    std::vector<double> expected = {3.0, 6.0, 9.0};
    EXPECT_EQ(result.size(), expected.size());
    for (size_t i = 0; i < result.size(); ++i) {
        EXPECT_DOUBLE_EQ(result[i], expected[i]);
    }
}

TEST_F(ReductionOpsTest, ReduceArrayMin) {
    ASSERT_MPI_PROCS(3);

    int rank = get_mpi_rank();
    std::vector<float> values = {static_cast<float>(rank + 10),
                                   static_cast<float>(rank + 20),
                                   static_cast<float>(rank + 30)};

    std::vector<float> result = ops->allreduce_array(values, ReductionOp::MIN);

    std::vector<float> expected = {10.0f, 20.0f, 30.0f};
    EXPECT_EQ(result.size(), expected.size());
    for (size_t i = 0; i < result.size(); ++i) {
        EXPECT_FLOAT_EQ(result[i], expected[i]);
    }
}

TEST_F(ReductionOpsTest, AllreduceArrayEmpty) {
    ASSERT_MPI_PROCS(2);

    std::vector<int> empty_array;
    std::vector<int> result = ops->allreduce_array(empty_array, ReductionOp::SUM);

    EXPECT_TRUE(result.empty());
}

// ============================================================================
// Broadcast Tests
// ============================================================================

TEST_F(ReductionOpsTest, BroadcastInt) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    int value = rank;  // Different value on each process

    ops->broadcast(value, 0);  // Broadcast from rank 0

    EXPECT_EQ(value, 0);  // All should have 0 now
}

TEST_F(ReductionOpsTest, BroadcastDouble) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    double value = rank * 3.14159;

    ops->broadcast(value, 2);  // Broadcast from rank 2

    EXPECT_DOUBLE_EQ(value, 2 * 3.14159);
}

TEST_F(ReductionOpsTest, BroadcastArray) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    std::vector<int> values;

    if (rank == 1) {
        values = {10, 20, 30, 40, 50};
    }

    ops->broadcast_array(values, 1);

    EXPECT_EQ(values.size(), 5);
    EXPECT_EQ(values[0], 10);
    EXPECT_EQ(values[1], 20);
    EXPECT_EQ(values[2], 30);
    EXPECT_EQ(values[3], 40);
    EXPECT_EQ(values[4], 50);
}

TEST_F(ReductionOpsTest, BroadcastArrayDouble) {
    ASSERT_MPI_PROCS(3);

    int rank = get_mpi_rank();
    std::vector<double> values;

    if (rank == 0) {
        values = {1.1, 2.2, 3.3, 4.4};
    }

    ops->broadcast_array(values, 0);

    EXPECT_EQ(values.size(), 4);
    EXPECT_DOUBLE_EQ(values[0], 1.1);
    EXPECT_DOUBLE_EQ(values[1], 2.2);
    EXPECT_DOUBLE_EQ(values[2], 3.3);
    EXPECT_DOUBLE_EQ(values[3], 4.4);
}

TEST_F(ReductionOpsTest, BroadcastEmptyArray) {
    ASSERT_MPI_PROCS(2);

    int rank = get_mpi_rank();
    std::vector<int> values;

    if (rank == 0) {
        // Keep it empty
    }

    ops->broadcast_array(values, 0);

    EXPECT_TRUE(values.empty());
}

// ============================================================================
// Barrier Tests
// ============================================================================

TEST_F(ReductionOpsTest, Barrier) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();

    // Artificially stagger processes
    if (rank == 0) {
        for (volatile int i = 0; i < 1000000; ++i) {}
    } else if (rank == 1) {
        for (volatile int i = 0; i < 500000; ++i) {}
    }

    ops->barrier();  // Synchronize

    // All processes should reach here
    EXPECT_TRUE(true);
}

// ============================================================================
// MaxLoc and MinLoc Tests
// ============================================================================

TEST_F(ReductionOpsTest, MaxLocInt) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    int value = rank * 10;  // 0, 10, 20, 30

    auto result = ops->maxloc(value);

    EXPECT_EQ(result.value, 30);
    EXPECT_EQ(result.rank, 3);
}

TEST_F(ReductionOpsTest, MinLocDouble) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    double value = (rank + 1) * 2.5;  // 2.5, 5.0, 7.5, 10.0

    auto result = ops->minloc(value);

    EXPECT_DOUBLE_EQ(result.value, 2.5);
    EXPECT_EQ(result.rank, 0);
}

TEST_F(ReductionOpsTest, MaxLocNegativeValues) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    int value = -rank;  // 0, -1, -2, -3

    auto result = ops->maxloc(value);

    EXPECT_EQ(result.value, 0);
    EXPECT_EQ(result.rank, 0);
}

TEST_F(ReductionOpsTest, MinLocLargeValues) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    int value = 1000 + rank;  // 1000, 1001, 1002, 1003

    auto result = ops->minloc(value);

    EXPECT_EQ(result.value, 1000);
    EXPECT_EQ(result.rank, 0);
}

// ============================================================================
// Different Types Tests
// ============================================================================

TEST_F(ReductionOpsTest, DifferentTypesInt) {
    ASSERT_MPI_PROCS(3);

    int rank = get_mpi_rank();
    int result = ops->allreduce(rank, ReductionOp::SUM);

    EXPECT_EQ(result, 0 + 1 + 2);
}

TEST_F(ReductionOpsTest, DifferentTypesFloat) {
    ASSERT_MPI_PROCS(3);

    int rank = get_mpi_rank();
    float value = static_cast<float>(rank) * 1.5f;
    float result = ops->allreduce(value, ReductionOp::SUM);

    EXPECT_FLOAT_EQ(result, 0.0f + 1.5f + 3.0f);
}

TEST_F(ReductionOpsTest, DifferentTypesDouble) {
    ASSERT_MPI_PROCS(3);

    int rank = get_mpi_rank();
    double value = rank * 2.5;
    double result = ops->allreduce(value, ReductionOp::SUM);

    EXPECT_DOUBLE_EQ(result, 0.0 + 2.5 + 5.0);
}

TEST_F(ReductionOpsTest, DifferentTypesUnsignedLong) {
    ASSERT_MPI_PROCS(3);

    int rank = get_mpi_rank();
    unsigned long value = static_cast<unsigned long>(rank * 100);
    unsigned long result = ops->allreduce(value, ReductionOp::SUM);

    EXPECT_EQ(result, 0 + 100 + 200);
}

// ============================================================================
// Custom Communicator Tests
// ============================================================================

TEST_F(ReductionOpsTest, CustomCommunicator) {
    ASSERT_MPI_PROCS(4);

    int rank = get_mpi_rank();
    int size = get_mpi_size();

    // Create a sub-communicator for even ranks
    MPI_Comm sub_comm;
    int color = (rank % 2 == 0) ? 0 : MPI_UNDEFINED;
    int key = rank;

    MPI_Comm_split(MPI_COMM_WORLD, color, key, &sub_comm);

    if (color == 0) {
        // Only even ranks participate
        ReductionOps sub_ops(sub_comm);

        int sub_rank;
        MPI_Comm_rank(sub_comm, &sub_rank);

        int result = sub_ops.allreduce(sub_rank, ReductionOp::SUM);

        // Ranks 0 and 2 are in sub-comm, with sub-ranks 0 and 1
        EXPECT_EQ(result, 0 + 1);
    }

    MPI_Comm_free(&sub_comm);
}

// ============================================================================
// MPIContext Integration Tests
// ============================================================================

TEST_F(ReductionOpsTest, MPIContextConstructor) {
    ASSERT_MPI_PROCS(3);

    int argc = 0;
    char** argv = nullptr;

    // Note: This test assumes MPI is already initialized
    // In a real scenario, you would use MPIContext first
    if (is_mpi_initialized()) {
        ReductionOps ctx_ops(MPI_COMM_WORLD);

        int rank = get_mpi_rank();
        int result = ctx_ops.allreduce(rank, ReductionOp::SUM);

        EXPECT_EQ(result, 0 + 1 + 2);
    }
}

// ============================================================================
// Main function for MPI-enabled tests
// ============================================================================

