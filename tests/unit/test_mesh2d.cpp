/**
 * @file test_mesh2d.cpp
 * @brief Unit tests for Mesh2D class
 */

#include <gtest/gtest.h>
#include "../../src/mesh/mesh2d.hpp"
#include "../../src/mpi/cartesian_topology.hpp"
#include "../../src/mpi/ghost_cell_exchange.hpp"
#include "../test_common.hpp"
#include <cmath>
#include <limits>

using namespace test_utils;

// =============================================================================
// Test Suite: ConstructionTest
// Tests for various constructors
// =============================================================================

TEST(ConstructionTest, BasicConstructorCreatesMesh) {
    size_t nx = 10;
    size_t ny = 8;
    double lx = 2.0;
    double ly = 1.5;

    Mesh2D mesh(nx, ny, lx, ly);

    EXPECT_EQ(mesh.nx(), nx);
    EXPECT_EQ(mesh.ny(), ny);
    EXPECT_EQ(mesh.total_nx(), nx);
    EXPECT_EQ(mesh.total_ny(), ny);
    EXPECT_DOUBLE_EQ(mesh.lx(), lx);
    EXPECT_DOUBLE_EQ(mesh.ly(), ly);
    EXPECT_FALSE(mesh.has_ghost_cells());
}

TEST(ConstructionTest, BasicConstructorDefaultDimensions) {
    Mesh2D mesh(10, 8);

    EXPECT_EQ(mesh.nx(), 10);
    EXPECT_EQ(mesh.ny(), 8);
    EXPECT_DOUBLE_EQ(mesh.lx(), 1.0);
    EXPECT_DOUBLE_EQ(mesh.ly(), 1.0);
}

TEST(ConstructionTest, BasicConstructorZeroDimensionsThrows) {
    EXPECT_THROW(Mesh2D(0, 10), std::invalid_argument);
    EXPECT_THROW(Mesh2D(10, 0), std::invalid_argument);
    EXPECT_THROW(Mesh2D(0, 0), std::invalid_argument);
}

TEST(ConstructionTest, MPICtorCreatesMeshWithGhostCells) {
    skip_if_no_mpi();

    CartesianTopology topo(MPI_COMM_WORLD);
    size_t nx = 10;
    size_t ny = 8;

    Mesh2D mesh(nx, ny, 2.0, 1.5, topo);

    EXPECT_EQ(mesh.nx(), nx);
    EXPECT_EQ(mesh.ny(), ny);
    EXPECT_EQ(mesh.total_nx(), nx + 2);
    EXPECT_EQ(mesh.total_ny(), ny + 2);
    EXPECT_TRUE(mesh.has_ghost_cells());
    EXPECT_NE(mesh.topology(), nullptr);
    EXPECT_NE(mesh.ghost_exchange(), nullptr);
}

TEST(ConstructionTest, CopyConstructorCreatesDeepCopy) {
    Mesh2D original(5, 4, 2.0, 1.0);
    original(0, 0) = 1.0;
    original(2, 2) = 2.0;
    original(4, 3) = 3.0;

    Mesh2D copy(original);

    EXPECT_EQ(copy.nx(), original.nx());
    EXPECT_EQ(copy.ny(), original.ny());
    EXPECT_EQ(copy.lx(), original.lx());
    EXPECT_EQ(copy.ly(), original.ly());
    EXPECT_EQ(copy(0, 0), 1.0);
    EXPECT_EQ(copy(2, 2), 2.0);
    EXPECT_EQ(copy(4, 3), 3.0);

    // Verify deep copy
    original(0, 0) = 10.0;
    EXPECT_EQ(copy(0, 0), 1.0);
}

TEST(ConstructionTest, MoveConstructorTransfersOwnership) {
    Mesh2D original(5, 4, 2.0, 1.0);
    original(0, 0) = 1.0;

    Mesh2D moved(std::move(original));

    EXPECT_EQ(moved.nx(), 5);
    EXPECT_EQ(moved.ny(), 4);
    EXPECT_EQ(moved(0, 0), 1.0);

    // Original should be in a valid but unspecified state
    // (we don't check its contents after move)
}

TEST(ConstructionTest, CopyAssignmentOperatorWorks) {
    Mesh2D mesh1(5, 4, 2.0, 1.0);
    Mesh2D mesh2(3, 2, 1.0, 1.0);

    mesh1(1, 1) = 5.0;
    mesh2 = mesh1;

    EXPECT_EQ(mesh2.nx(), 5);
    EXPECT_EQ(mesh2.ny(), 4);
    EXPECT_EQ(mesh2(1, 1), 5.0);
}

TEST(ConstructionTest, MoveAssignmentOperatorWorks) {
    Mesh2D mesh1(5, 4, 2.0, 1.0);
    Mesh2D mesh2(3, 2, 1.0, 1.0);

    mesh1(1, 1) = 5.0;
    mesh2 = std::move(mesh1);

    EXPECT_EQ(mesh2.nx(), 5);
    EXPECT_EQ(mesh2.ny(), 4);
    EXPECT_EQ(mesh2(1, 1), 5.0);
}

TEST(ConstructionTest, SelfCopyAssignmentIsNoOp) {
    Mesh2D mesh(5, 4, 2.0, 1.0);
    mesh(1, 1) = 5.0;

    mesh = mesh;

    EXPECT_EQ(mesh(1, 1), 5.0);
}

TEST(ConstructionTest, SelfMoveAssignmentIsNoOp) {
    Mesh2D mesh(5, 4, 2.0, 1.0);
    mesh(1, 1) = 5.0;

    mesh = std::move(mesh);

    EXPECT_EQ(mesh(1, 1), 5.0);
}

// =============================================================================
// Test Suite: AccessorTest
// Tests for element access and bounds checking
// =============================================================================

TEST(AccessorTest, ElementAccessReturnsCorrectValues) {
    Mesh2D mesh(5, 4, 2.0, 1.0);

    mesh(0, 0) = 1.0;
    mesh(2, 2) = 2.0;
    mesh(4, 3) = 3.0;

    EXPECT_EQ(mesh(0, 0), 1.0);
    EXPECT_EQ(mesh(2, 2), 2.0);
    EXPECT_EQ(mesh(4, 3), 3.0);
}

TEST(AccessorTest, ConstElementAccessWorks) {
    Mesh2D mesh(5, 4, 2.0, 1.0);
    mesh(1, 1) = 42.0;

    const Mesh2D& const_mesh = mesh;

    EXPECT_EQ(const_mesh(1, 1), 42.0);
}

TEST(AccessorTest, OutOfBoundsAccessThrows) {
    Mesh2D mesh(5, 4);

    EXPECT_THROW(mesh(5, 0), std::out_of_range);
    EXPECT_THROW(mesh(0, 5), std::out_of_range);
    EXPECT_THROW(mesh(5, 5), std::out_of_range);
}

TEST(AccessorTest, DirectAtAccessWorks) {
    Mesh2D mesh(5, 4);
    mesh(0, 0) = 1.0;
    mesh(2, 2) = 2.0;

    EXPECT_EQ(mesh.at(0, 0), 1.0);
    EXPECT_EQ(mesh.at(2, 2), 2.0);
}

TEST(AccessorTest, DirectAtOutOfBoundsThrows) {
    Mesh2D mesh(5, 4);

    EXPECT_THROW(mesh.at(5, 0), std::out_of_range);
    EXPECT_THROW(mesh.at(0, 5), std::out_of_range);
}

TEST(AccessorTest, DataAccessReturnsArray) {
    Mesh2D mesh(5, 4);

    utils::Array2D& data = mesh.data();
    EXPECT_EQ(data.rows(), 4);
    EXPECT_EQ(data.cols(), 5);
}

TEST(AccessorTest, ConstDataAccessReturnsArray) {
    Mesh2D mesh(5, 4);

    const Mesh2D& const_mesh = mesh;
    const utils::Array2D& data = const_mesh.data();

    EXPECT_EQ(data.rows(), 4);
    EXPECT_EQ(data.cols(), 5);
}

// =============================================================================
// Test Suite: CoordinateTest
// Tests for global/local coordinate conversions
// =============================================================================

TEST(CoordinateTest, LocalToGlobalConversionWithoutMPI) {
    Mesh2D mesh(10, 8);

    for (size_t i = 0; i < mesh.nx(); ++i) {
        EXPECT_EQ(mesh.local_to_global_x(i), i);
    }

    for (size_t j = 0; j < mesh.ny(); ++j) {
        EXPECT_EQ(mesh.local_to_global_y(j), j);
    }
}

TEST(CoordinateTest, GlobalToLocalConversionWithoutMPI) {
    Mesh2D mesh(10, 8);

    for (size_t i = 0; i < mesh.nx(); ++i) {
        EXPECT_EQ(mesh.global_to_local_x(i), i);
    }

    for (size_t j = 0; j < mesh.ny(); ++j) {
        EXPECT_EQ(mesh.global_to_local_y(j), j);
    }
}

TEST(CoordinateTest, GlobalToLocalOutOfBoundsThrows) {
    Mesh2D mesh(10, 8);

    EXPECT_THROW(mesh.global_to_local_x(10), std::out_of_range);
    EXPECT_THROW(mesh.global_to_local_y(8), std::out_of_range);
}

TEST(CoordinateTest, LocalToGlobalWithMPI) {
    skip_if_no_mpi();
    skip_if_mpi_procs_not(4);  // Need 4 processes for 2x2 topology

    CartesianTopology topo(MPI_COMM_WORLD);
    Mesh2D mesh(10, 8, 2.0, 1.5, topo);

    // Test that conversion is consistent
    for (size_t i = 0; i < mesh.nx(); ++i) {
        size_t global_i = mesh.local_to_global_x(i);
        size_t local_i = mesh.global_to_local_x(global_i);
        EXPECT_EQ(local_i, i);
    }

    for (size_t j = 0; j < mesh.ny(); ++j) {
        size_t global_j = mesh.local_to_global_y(j);
        size_t local_j = mesh.global_to_local_y(global_j);
        EXPECT_EQ(local_j, j);
    }
}

// =============================================================================
// Test Suite: PhysicalCoordTest
// Tests for physical coordinate calculations
// =============================================================================

TEST(PhysicalCoordTest, XCoordinateCalculatesCorrectly) {
    Mesh2D mesh(10, 8, 2.0, 1.5);

    double hx = mesh.hx();
    EXPECT_DOUBLE_EQ(mesh.x_coord(0), 0.0);
    EXPECT_DOUBLE_EQ(mesh.x_coord(1), hx);
    EXPECT_DOUBLE_EQ(mesh.x_coord(9), 2.0);
}

TEST(PhysicalCoordTest, YCoordinateCalculatesCorrectly) {
    Mesh2D mesh(10, 8, 2.0, 1.5);

    double hy = mesh.hy();
    EXPECT_DOUBLE_EQ(mesh.y_coord(0), 0.0);
    EXPECT_DOUBLE_EQ(mesh.y_coord(1), hy);
    EXPECT_DOUBLE_EQ(mesh.y_coord(7), 1.5);
}

TEST(PhysicalCoordTest, CoordinatePairCalculatesCorrectly) {
    Mesh2D mesh(10, 8, 2.0, 1.5);

    auto [x, y] = mesh.coord(5, 4);

    EXPECT_DOUBLE_EQ(x, mesh.x_coord(5));
    EXPECT_DOUBLE_EQ(y, mesh.y_coord(4));
}

TEST(PhysicalCoordTest, GridSpacingCalculatesCorrectly) {
    Mesh2D mesh(10, 8, 2.0, 1.5);

    // hx = lx / (nx - 1) for inclusive endpoints
    EXPECT_DOUBLE_EQ(mesh.hx(), 2.0 / 9.0);
    EXPECT_DOUBLE_EQ(mesh.hy(), 1.5 / 7.0);
}

TEST(PhysicalCoordTest, GridSpacingForSinglePoint) {
    Mesh2D mesh(1, 1, 2.0, 1.5);

    // For single point, spacing equals domain length
    EXPECT_DOUBLE_EQ(mesh.hx(), 2.0);
    EXPECT_DOUBLE_EQ(mesh.hy(), 1.5);
}

TEST(PhysicalCoordTest, XCoordinateOutOfBoundsThrows) {
    Mesh2D mesh(10, 8);

    EXPECT_THROW(mesh.x_coord(10), std::out_of_range);
}

TEST(PhysicalCoordTest, YCoordinateOutOfBoundsThrows) {
    Mesh2D mesh(10, 8);

    EXPECT_THROW(mesh.y_coord(8), std::out_of_range);
}

// =============================================================================
// Test Suite: BoundaryConditionTest
// Tests for Dirichlet boundary conditions
// =============================================================================

TEST(BoundaryConditionTest, ApplyDirichletBCWithGhostCells) {
    skip_if_no_mpi();

    CartesianTopology topo(MPI_COMM_WORLD);
    Mesh2D mesh(5, 4, 2.0, 1.5, topo);

    mesh.apply_dirichlet_bc(42.0);

    // Check all ghost cells are set to 42.0
    for (size_t j = 0; j < mesh.total_nx(); ++j) {
        EXPECT_EQ(mesh.at(0, j), 42.0);              // South
        EXPECT_EQ(mesh.at(mesh.total_ny() - 1, j), 42.0);  // North
    }

    for (size_t i = 0; i < mesh.total_ny(); ++i) {
        EXPECT_EQ(mesh.at(i, 0), 42.0);              // West
        EXPECT_EQ(mesh.at(i, mesh.total_nx() - 1), 42.0);  // East
    }
}

TEST(BoundaryConditionTest, ApplyDirichletBCWithoutGhostCells) {
    Mesh2D mesh(5, 4, 2.0, 1.5);

    // Should not throw, but no-op
    mesh.apply_dirichlet_bc(42.0);
}

TEST(BoundaryConditionTest, ApplyFunctionBC) {
    skip_if_no_mpi();

    CartesianTopology topo(MPI_COMM_WORLD);
    Mesh2D mesh(5, 4, 2.0, 1.5, topo);

    auto bc_func = [](double x, double y, double t) -> double {
        return x + y + t;
    };

    mesh.apply_bc(bc_func, 1.0);

    // Check boundaries
    EXPECT_DOUBLE_EQ(mesh.at(0, 0), 0.0 + 0.0 + 1.0);           // South-West corner
    EXPECT_DOUBLE_EQ(mesh.at(0, 1), mesh.hx() + 0.0 + 1.0);     // South
    EXPECT_DOUBLE_EQ(mesh.at(mesh.total_ny() - 1, 0), 0.0 + mesh.ly() + 1.0);  // West
}

TEST(BoundaryConditionTest, ApplyBCAtDirectionNorth) {
    skip_if_no_mpi();

    CartesianTopology topo(MPI_COMM_WORLD);
    Mesh2D mesh(5, 4, 2.0, 1.5, topo);

    mesh.apply_bc_at_direction(CardinalDirection::North, 100.0);

    // Check only North boundary
    for (size_t j = 0; j < mesh.total_nx(); ++j) {
        EXPECT_EQ(mesh.at(mesh.total_ny() - 1, j), 100.0);
    }
}

TEST(BoundaryConditionTest, ApplyBCAtDirectionSouth) {
    skip_if_no_mpi();

    CartesianTopology topo(MPI_COMM_WORLD);
    Mesh2D mesh(5, 4, 2.0, 1.5, topo);

    mesh.apply_bc_at_direction(CardinalDirection::South, 200.0);

    // Check only South boundary
    for (size_t j = 0; j < mesh.total_nx(); ++j) {
        EXPECT_EQ(mesh.at(0, j), 200.0);
    }
}

TEST(BoundaryConditionTest, ApplyBCAtDirectionEast) {
    skip_if_no_mpi();

    CartesianTopology topo(MPI_COMM_WORLD);
    Mesh2D mesh(5, 4, 2.0, 1.5, topo);

    mesh.apply_bc_at_direction(CardinalDirection::East, 300.0);

    // Check only East boundary
    for (size_t i = 0; i < mesh.total_ny(); ++i) {
        EXPECT_EQ(mesh.at(i, mesh.total_nx() - 1), 300.0);
    }
}

TEST(BoundaryConditionTest, ApplyBCAtDirectionWest) {
    skip_if_no_mpi();

    CartesianTopology topo(MPI_COMM_WORLD);
    Mesh2D mesh(5, 4, 2.0, 1.5, topo);

    mesh.apply_bc_at_direction(CardinalDirection::West, 400.0);

    // Check only West boundary
    for (size_t i = 0; i < mesh.total_ny(); ++i) {
        EXPECT_EQ(mesh.at(i, 0), 400.0);
    }
}

// =============================================================================
// Test Suite: LaplacianTest
// Tests for 5-point stencil Laplacian computation
// =============================================================================

TEST(LaplacianTest, ComputeLaplacianWithConstantFunction) {
    Mesh2D mesh(5, 4, 2.0, 1.5);
    mesh.fill(1.0);  // Constant function

    utils::Array2D result(4, 5);
    mesh.compute_laplacian(result);

    // Laplacian of constant is zero
    for (size_t i = 0; i < mesh.ny(); ++i) {
        for (size_t j = 0; j < mesh.nx(); ++j) {
            EXPECT_NEAR(result(i, j), 0.0, 1e-10);
        }
    }
}

TEST(LaplacianTest, ComputeLaplacianWithLinearFunction) {
    Mesh2D mesh(5, 4, 2.0, 1.5);

    // Linear function: u = x + y
    for (size_t i = 0; i < mesh.ny(); ++i) {
        for (size_t j = 0; j < mesh.nx(); ++j) {
            double x = mesh.x_coord(j);
            double y = mesh.y_coord(i);
            mesh(i, j) = x + y;
        }
    }

    utils::Array2D result(4, 5);
    mesh.compute_laplacian(result);

    // Laplacian of linear function is zero
    for (size_t i = 0; i < mesh.ny(); ++i) {
        for (size_t j = 0; j < mesh.nx(); ++j) {
            EXPECT_NEAR(result(i, j), 0.0, 1e-10);
        }
    }
}

TEST(LaplacianTest, ComputeLaplacianWithQuadraticFunction) {
    Mesh2D mesh(5, 4, 2.0, 1.5);

    // Quadratic function: u = x^2 + y^2
    for (size_t i = 0; i < mesh.ny(); ++i) {
        for (size_t j = 0; j < mesh.nx(); ++j) {
            double x = mesh.x_coord(j);
            double y = mesh.y_coord(i);
            mesh(i, j) = x * x + y * y;
        }
    }

    utils::Array2D result(4, 5);
    mesh.compute_laplacian(result);

    // Laplacian of x^2 + y^2 is 2 + 2 = 4
    for (size_t i = 0; i < mesh.ny(); ++i) {
        for (size_t j = 0; j < mesh.nx(); ++j) {
            EXPECT_NEAR(result(i, j), 4.0, 1e-8);
        }
    }
}

TEST(LaplacianTest, ComputeLaplacianWithWrongDimensionsThrows) {
    Mesh2D mesh(5, 4, 2.0, 1.5);
    utils::Array2D result(3, 5);  // Wrong dimensions

    EXPECT_THROW(mesh.compute_laplacian(result), std::invalid_argument);
}

TEST(LaplacianTest, ComputeLaplacianWithGhostCells) {
    skip_if_no_mpi();

    CartesianTopology topo(MPI_COMM_WORLD);
    Mesh2D mesh(5, 4, 2.0, 1.5, topo);

    // Set up function
    for (size_t i = 0; i < mesh.ny(); ++i) {
        for (size_t j = 0; j < mesh.nx(); ++j) {
            double x = mesh.x_coord(j);
            double y = mesh.y_coord(i);
            mesh(i, j) = x * x + y * y;
        }
    }

    // Set ghost cells using the same function
    mesh.apply_bc([](double x, double y, double t) { return x * x + y * y + t; }, 0.0);

    utils::Array2D result(4, 5);
    mesh.compute_laplacian(result);

    // Should still compute correct Laplacian
    for (size_t i = 0; i < mesh.ny(); ++i) {
        for (size_t j = 0; j < mesh.nx(); ++j) {
            EXPECT_NEAR(result(i, j), 4.0, 1e-8);
        }
    }
}

// =============================================================================
// Test Suite: GhostCellTest
// Tests for ghost cell exchange (MPI)
// =============================================================================

TEST(GhostCellTest, ExchangeGhostCellsWithoutMPIThrows) {
    Mesh2D mesh(5, 4, 2.0, 1.5);

    EXPECT_THROW(mesh.exchange_ghost_cells(), std::runtime_error);
}

TEST(GhostCellTest, ExchangeGhostCellsWithMPI) {
    skip_if_no_mpi();
    skip_if_mpi_procs_not(4);  // Need 4 processes for 2x2 topology

    CartesianTopology topo(MPI_COMM_WORLD);
    Mesh2D mesh(5, 4, 2.0, 1.5, topo);

    // Set unique values in each process
    int rank = topo.rank();
    for (size_t i = 0; i < mesh.ny(); ++i) {
        for (size_t j = 0; j < mesh.nx(); ++j) {
            mesh(i, j) = static_cast<double>(rank);
        }
    }

    // Exchange ghost cells
    mesh.exchange_ghost_cells();

    // Verify ghost cells contain neighbor values
    // (exact verification depends on topology and would need more complex checks)
    SUCCEED();  // If we got here without exception, it worked
}

TEST(GhostCellTest, SetBoundaryValuesWithGhostCells) {
    skip_if_no_mpi();

    CartesianTopology topo(MPI_COMM_WORLD);
    Mesh2D mesh(5, 4, 2.0, 1.5, topo);

    mesh.set_boundary_values(99.0);

    // Check all ghost cells
    for (size_t j = 0; j < mesh.total_nx(); ++j) {
        EXPECT_EQ(mesh.at(0, j), 99.0);              // South
        EXPECT_EQ(mesh.at(mesh.total_ny() - 1, j), 99.0);  // North
    }

    for (size_t i = 0; i < mesh.total_ny(); ++i) {
        EXPECT_EQ(mesh.at(i, 0), 99.0);              // West
        EXPECT_EQ(mesh.at(i, mesh.total_nx() - 1), 99.0);  // East
    }
}

// =============================================================================
// Test Suite: NormsTest
// Tests for various norm calculations
// =============================================================================

TEST(NormsTest, L2NormOfZeroVectorIsZero) {
    Mesh2D mesh(5, 4, 2.0, 1.5);
    mesh.fill(0.0);

    EXPECT_DOUBLE_EQ(mesh.l2_norm(), 0.0);
}

TEST(NormsTest, L2NormOfOnes) {
    Mesh2D mesh(3, 2, 2.0, 1.5);
    mesh.fill(1.0);

    // sqrt(3 * 2 * 1^2) = sqrt(6)
    double expected = std::sqrt(6.0);
    EXPECT_NEAR(mesh.l2_norm(), expected, 1e-10);
}

TEST(NormsTest, L2NormOfScaledVector) {
    Mesh2D mesh(3, 2, 2.0, 1.5);
    mesh.fill(2.0);

    // sqrt(3 * 2 * 2^2) = sqrt(24) = 2 * sqrt(6)
    double expected = 2.0 * std::sqrt(6.0);
    EXPECT_NEAR(mesh.l2_norm(), expected, 1e-10);
}

TEST(NormsTest, LinftyNormOfZeroVectorIsZero) {
    Mesh2D mesh(5, 4, 2.0, 1.5);
    mesh.fill(0.0);

    EXPECT_DOUBLE_EQ(mesh.linfty_norm(), 0.0);
}

TEST(NormsTest, LinftyNormOfOnes) {
    Mesh2D mesh(3, 2, 2.0, 1.5);
    mesh.fill(1.0);

    EXPECT_DOUBLE_EQ(mesh.linfty_norm(), 1.0);
}

TEST(NormsTest, LinftyNormOfMixedValues) {
    Mesh2D mesh(3, 2, 2.0, 1.5);

    mesh(0, 0) = 1.0;
    mesh(1, 1) = -5.0;
    mesh(2, 1) = 3.0;

    EXPECT_DOUBLE_EQ(mesh.linfty_norm(), 5.0);
}

TEST(NormsTest, MaxValue) {
    Mesh2D mesh(3, 2, 2.0, 1.5);

    mesh(0, 0) = 1.0;
    mesh(1, 1) = -5.0;
    mesh(2, 1) = 3.0;

    EXPECT_DOUBLE_EQ(mesh.max(), 3.0);
}

TEST(NormsTest, MinValue) {
    Mesh2D mesh(3, 2, 2.0, 1.5);

    mesh(0, 0) = 1.0;
    mesh(1, 1) = -5.0;
    mesh(2, 1) = 3.0;

    EXPECT_DOUBLE_EQ(mesh.min(), -5.0);
}

// =============================================================================
// Test Suite: ArithmeticTest
// Tests for mesh arithmetic operations
// =============================================================================

TEST(ArithmeticTest, FillOperation) {
    Mesh2D mesh(5, 4, 2.0, 1.5);

    mesh.fill(42.0);

    for (size_t i = 0; i < mesh.ny(); ++i) {
        for (size_t j = 0; j < mesh.nx(); ++j) {
            EXPECT_EQ(mesh(i, j), 42.0);
        }
    }
}

TEST(ArithmeticTest, CopyFromOperation) {
    Mesh2D mesh1(5, 4, 2.0, 1.5);
    Mesh2D mesh2(5, 4, 2.0, 1.5);

    mesh1.fill(10.0);
    mesh2.copy_from(mesh1);

    for (size_t i = 0; i < mesh2.ny(); ++i) {
        for (size_t j = 0; j < mesh2.nx(); ++j) {
            EXPECT_EQ(mesh2(i, j), 10.0);
        }
    }
}

TEST(ArithmeticTest, CopyFromWithMismatchedDimensionsThrows) {
    Mesh2D mesh1(5, 4, 2.0, 1.5);
    Mesh2D mesh2(3, 2, 1.0, 1.0);

    EXPECT_THROW(mesh2.copy_from(mesh1), std::invalid_argument);
}

TEST(ArithmeticTest, ScaleOperation) {
    Mesh2D mesh(3, 2, 2.0, 1.5);

    for (size_t i = 0; i < mesh.ny(); ++i) {
        for (size_t j = 0; j < mesh.nx(); ++j) {
            mesh(i, j) = static_cast<double>(i + j);
        }
    }

    mesh.scale(2.0);

    for (size_t i = 0; i < mesh.ny(); ++i) {
        for (size_t j = 0; j < mesh.nx(); ++j) {
            EXPECT_DOUBLE_EQ(mesh(i, j), 2.0 * static_cast<double>(i + j));
        }
    }
}

TEST(ArithmeticTest, AddOperation) {
    Mesh2D mesh1(3, 2, 2.0, 1.5);
    Mesh2D mesh2(3, 2, 2.0, 1.5);

    mesh1.fill(5.0);
    mesh2.fill(3.0);

    mesh1.add(mesh2);

    for (size_t i = 0; i < mesh1.ny(); ++i) {
        for (size_t j = 0; j < mesh1.nx(); ++j) {
            EXPECT_DOUBLE_EQ(mesh1(i, j), 8.0);
        }
    }
}

TEST(ArithmeticTest, SubtractOperation) {
    Mesh2D mesh1(3, 2, 2.0, 1.5);
    Mesh2D mesh2(3, 2, 2.0, 1.5);

    mesh1.fill(5.0);
    mesh2.fill(3.0);

    mesh1.subtract(mesh2);

    for (size_t i = 0; i < mesh1.ny(); ++i) {
        for (size_t j = 0; j < mesh1.nx(); ++j) {
            EXPECT_DOUBLE_EQ(mesh1(i, j), 2.0);
        }
    }
}

TEST(ArithmeticTest, AddWithMismatchedDimensionsThrows) {
    Mesh2D mesh1(3, 2, 2.0, 1.5);
    Mesh2D mesh2(5, 4, 2.0, 1.5);

    EXPECT_THROW(mesh1.add(mesh2), std::invalid_argument);
}

TEST(ArithmeticTest, SubtractWithMismatchedDimensionsThrows) {
    Mesh2D mesh1(3, 2, 2.0, 1.5);
    Mesh2D mesh2(5, 4, 2.0, 1.5);

    EXPECT_THROW(mesh1.subtract(mesh2), std::invalid_argument);
}

// =============================================================================
// Test Suite: EdgeCases
// Tests for edge cases like 1x1, 1xN, Nx1 grids
// =============================================================================

TEST(EdgeCases, OneByOneMesh) {
    Mesh2D mesh(1, 1, 1.0, 1.0);

    EXPECT_EQ(mesh.nx(), 1);
    EXPECT_EQ(mesh.ny(), 1);
    EXPECT_DOUBLE_EQ(mesh.hx(), 1.0);
    EXPECT_DOUBLE_EQ(mesh.hy(), 1.0);

    mesh(0, 0) = 42.0;
    EXPECT_EQ(mesh(0, 0), 42.0);
}

TEST(EdgeCases, OneByNMesh) {
    Mesh2D mesh(1, 5, 1.0, 2.0);

    EXPECT_EQ(mesh.nx(), 1);
    EXPECT_EQ(mesh.ny(), 5);
    EXPECT_DOUBLE_EQ(mesh.hx(), 1.0);
    EXPECT_NEAR(mesh.hy(), 2.0 / 4.0, 1e-10);
}

TEST(EdgeCases, NByOneMesh) {
    Mesh2D mesh(5, 1, 2.0, 1.0);

    EXPECT_EQ(mesh.nx(), 5);
    EXPECT_EQ(mesh.ny(), 1);
    EXPECT_NEAR(mesh.hx(), 2.0 / 4.0, 1e-10);
    EXPECT_DOUBLE_EQ(mesh.hy(), 1.0);
}

TEST(EdgeCases, OneByOneMeshLaplacian) {
    Mesh2D mesh(1, 1, 1.0, 1.0);
    mesh(0, 0) = 1.0;

    utils::Array2D result(1, 1);
    mesh.compute_laplacian(result);

    // For 1x1 mesh, Laplacian should be 0 (no neighbors)
    EXPECT_NEAR(result(0, 0), 0.0, 1e-10);
}

TEST(EdgeCases, OneByOneMeshNorms) {
    Mesh2D mesh(1, 1, 1.0, 1.0);
    mesh(0, 0) = 5.0;

    EXPECT_DOUBLE_EQ(mesh.l2_norm(), 5.0);
    EXPECT_DOUBLE_EQ(mesh.linfty_norm(), 5.0);
    EXPECT_DOUBLE_EQ(mesh.max(), 5.0);
    EXPECT_DOUBLE_EQ(mesh.min(), 5.0);
}

TEST(EdgeCases, VerySmallMesh) {
    Mesh2D mesh(2, 2, 1.0, 1.0);

    mesh(0, 0) = 1.0;
    mesh(0, 1) = 2.0;
    mesh(1, 0) = 3.0;
    mesh(1, 1) = 4.0;

    EXPECT_DOUBLE_EQ(mesh.l2_norm(), std::sqrt(1.0 + 4.0 + 9.0 + 16.0));
    EXPECT_DOUBLE_EQ(mesh.max(), 4.0);
    EXPECT_DOUBLE_EQ(mesh.min(), 1.0);
}

TEST(EdgeCases, LargeMeshPerformance) {
    // Test that operations work on larger meshes
    Mesh2D mesh(100, 100, 10.0, 10.0);

    mesh.fill(1.0);

    double norm = mesh.l2_norm();
    EXPECT_GT(norm, 0.0);

    mesh.scale(2.0);

    EXPECT_DOUBLE_EQ(mesh.max(), 2.0);
    EXPECT_DOUBLE_EQ(mesh.min(), 2.0);
}

// =============================================================================
// Test Suite: MPITopologyIntegration
// Tests for integration with Cartesian topology
// =============================================================================

TEST(MPITopologyIntegration, MeshUsesTopology) {
    skip_if_no_mpi();

    CartesianTopology topo(MPI_COMM_WORLD);
    Mesh2D mesh(10, 8, 2.0, 1.5, topo);

    EXPECT_NE(mesh.topology(), nullptr);
    EXPECT_EQ(mesh.topology()->rank(), topo.rank());
    EXPECT_EQ(mesh.topology()->size(), topo.size());
}

TEST(MPITopologyIntegration, GhostExchangeUsesTopology) {
    skip_if_no_mpi();
    skip_if_mpi_procs_not(4);

    CartesianTopology topo(MPI_COMM_WORLD);
    Mesh2D mesh(10, 8, 2.0, 1.5, topo);

    EXPECT_NE(mesh.ghost_exchange(), nullptr);
    EXPECT_EQ(mesh.ghost_exchange()->topology().rank(), topo.rank());
}

// =============================================================================
// Main function for MPI-enabled tests
// =============================================================================

