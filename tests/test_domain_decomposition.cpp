/**
 * @file test_domain_decomposition.cpp
 * @brief Unit tests for DomainDecomposition class
 */

#include <gtest/gtest.h>
#include "mesh/domain_decomposition.hpp"

// Mock CartesianTopology for testing
class MockCartesianTopology {
public:
    MockCartesianTopology(int rank, int size, int dim_x, int dim_y)
        : rank_(rank), size_(size), dims_({dim_x, dim_y}),
          coords_({rank % dim_x, rank / dim_x}) {}

    int rank() const { return rank_; }
    int size() const { return size_; }
    MPI_Comm communicator() const { return MPI_COMM_WORLD; }
    const std::vector<int>& dims() const { return dims_; }
    const std::vector<int>& coords() const { return coords_; }

    bool is_on_boundary(Direction dir, Shift shift) const {
        switch (dir) {
            case Direction::X:
                if (shift == Shift::Backward) return coords_[0] == 0;
                else return coords_[0] == dims_[0] - 1;
            case Direction::Y:
                if (shift == Shift::Backward) return coords_[1] == 0;
                else return coords_[1] == dims_[1] - 1;
        }
        return false;
    }

private:
    int rank_;
    int size_;
    std::vector<int> dims_;
    std::vector<int> coords_;
};

// Test compute_optimal_dimensions
TEST(DomainDecompositionTest, ComputeOptimalDimensions) {
    // Test with 1 process
    auto dims1 = DomainDecomposition::compute_optimal_dimensions(1, 100, 100);
    EXPECT_EQ(dims1.first, 1);
    EXPECT_EQ(dims1.second, 1);

    // Test with 4 processes (perfect square)
    auto dims4 = DomainDecomposition::compute_optimal_dimensions(4, 100, 100);
    EXPECT_EQ(dims4.first, 2);
    EXPECT_EQ(dims4.second, 2);

    // Test with 6 processes (2x3)
    auto dims6 = DomainDecomposition::compute_optimal_dimensions(6, 100, 100);
    EXPECT_EQ(dims6.first, 2);
    EXPECT_EQ(dims6.second, 3);

    // Test with 12 processes (3x4)
    auto dims12 = DomainDecomposition::compute_optimal_dimensions(12, 100, 100);
    EXPECT_EQ(dims12.first, 3);
    EXPECT_EQ(dims12.second, 4);
}

// Test Subdomain coordinate transformation
TEST(SubdomainTest, CoordinateTransformation) {
    Subdomain subdomain = {10, 20, 15, 25, 100, 100};

    // Test valid coordinate transformation
    size_t local_i = subdomain.global_to_local_i(15);
    EXPECT_EQ(local_i, 5);

    size_t local_j = subdomain.global_to_local_j(25);
    EXPECT_EQ(local_j, 5);

    // Test out of bounds throws exception
    EXPECT_THROW(subdomain.global_to_local_i(5), std::out_of_range);
    EXPECT_THROW(subdomain.global_to_local_i(30), std::out_of_range);
    EXPECT_THROW(subdomain.global_to_local_j(10), std::out_of_range);
    EXPECT_THROW(subdomain.global_to_local_j(50), std::out_of_range);
}

// Test manual decomposition validation
TEST(DomainDecompositionTest, ManualDecompositionValidation) {
    MockCartesianTopology topology(0, 4, 2, 2);

    // Valid non-overlapping decomposition
    std::vector<Subdomain> subdomains = {
        {0, 0, 50, 50, 100, 100},
        {50, 0, 50, 50, 100, 100},
        {0, 50, 50, 50, 100, 100},
        {50, 50, 50, 50, 100, 100}
    };

    EXPECT_NO_THROW({
        DomainDecomposition decomp(subdomains, topology);
    });

    // Invalid: overlapping subdomains
    std::vector<Subdomain> overlapping = {
        {0, 0, 60, 50, 100, 100},  // Overlaps with next
        {50, 0, 50, 50, 100, 100},
        {0, 50, 50, 50, 100, 100},
        {50, 50, 50, 50, 100, 100}
    };

    EXPECT_THROW({
        DomainDecomposition decomp(overlapping, topology);
    }, std::invalid_argument);

    // Invalid: incomplete coverage
    std::vector<Subdomain> incomplete = {
        {0, 0, 40, 50, 100, 100},
        {50, 0, 50, 50, 100, 100},
        {0, 50, 50, 50, 100, 100},
        {50, 50, 50, 50, 100, 100}
    };

    EXPECT_THROW({
        DomainDecomposition decomp(incomplete, topology);
    }, std::invalid_argument);
}
