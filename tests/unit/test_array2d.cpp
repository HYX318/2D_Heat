/**
 * @file test_array2d.cpp
 * @brief Unit tests for Array2D class using Google Test framework
 */

#include <gtest/gtest.h>
#include "utils/array2d.hpp"
#include <cmath>
#include <limits>

using utils::Array2D;

/**
 * @test Construction - Test default constructor
 */
TEST(Array2DTest, Construction) {
    // Test basic construction
    Array2D arr(3, 4);
    EXPECT_EQ(arr.rows(), 3);
    EXPECT_EQ(arr.cols(), 4);
    EXPECT_EQ(arr.size(), 12);
    EXPECT_NE(arr.data(), nullptr);

    // Test square array
    Array2D square(5, 5);
    EXPECT_EQ(square.rows(), 5);
    EXPECT_EQ(square.cols(), 5);
    EXPECT_EQ(square.size(), 25);

    // Test invalid dimensions
    EXPECT_THROW(Array2D(0, 5), std::invalid_argument);
    EXPECT_THROW(Array2D(5, 0), std::invalid_argument);
}

/**
 * @test FillConstruction - Test fill constructor
 */
TEST(Array2DTest, FillConstruction) {
    Array2D arr(3, 4, 2.5);

    EXPECT_EQ(arr.rows(), 3);
    EXPECT_EQ(arr.cols(), 4);
    EXPECT_EQ(arr.size(), 12);

    // Verify all elements are filled with the specified value
    for (size_t i = 0; i < arr.rows(); ++i) {
        for (size_t j = 0; j < arr.cols(); ++j) {
            EXPECT_DOUBLE_EQ(arr(i, j), 2.5);
        }
    }

    // Test with zero value
    Array2D zeros(2, 3, 0.0);
    for (size_t i = 0; i < zeros.rows(); ++i) {
        for (size_t j = 0; j < zeros.cols(); ++j) {
            EXPECT_DOUBLE_EQ(zeros(i, j), 0.0);
        }
    }

    // Test with negative value
    Array2D negatives(2, 2, -3.14);
    for (size_t i = 0; i < negatives.rows(); ++i) {
        for (size_t j = 0; j < negatives.cols(); ++j) {
            EXPECT_DOUBLE_EQ(negatives(i, j), -3.14);
        }
    }
}

/**
 * @test CopyConstruction - Test deep copy
 */
TEST(Array2DTest, CopyConstruction) {
    // Create original array and fill with values
    Array2D original(3, 3);
    double counter = 1.0;
    for (size_t i = 0; i < original.rows(); ++i) {
        for (size_t j = 0; j < original.cols(); ++j) {
            original(i, j) = counter++;
        }
    }

    // Create copy
    Array2D copy(original);

    // Verify dimensions
    EXPECT_EQ(copy.rows(), original.rows());
    EXPECT_EQ(copy.cols(), original.cols());
    EXPECT_EQ(copy.size(), original.size());

    // Verify deep copy (same values, different memory)
    for (size_t i = 0; i < original.rows(); ++i) {
        for (size_t j = 0; j < original.cols(); ++j) {
            EXPECT_DOUBLE_EQ(copy(i, j), original(i, j));
        }
    }
    EXPECT_NE(copy.data(), original.data());

    // Modify original and verify copy is independent
    original(0, 0) = 999.0;
    EXPECT_DOUBLE_EQ(copy(0, 0), 1.0);
}

/**
 * @test MoveConstruction - Test move semantics
 */
TEST(Array2DTest, MoveConstruction) {
    // Create original array
    Array2D original(3, 4, 5.0);
    double* original_data = original.data();
    size_t original_rows = original.rows();
    size_t original_cols = original.cols();

    // Move construct
    Array2D moved(std::move(original));

    // Verify moved array has correct properties
    EXPECT_EQ(moved.rows(), original_rows);
    EXPECT_EQ(moved.cols(), original_cols);
    EXPECT_EQ(moved.size(), original_rows * original_cols);
    EXPECT_EQ(moved.data(), original_data);

    // Verify original array is in valid but empty state
    EXPECT_EQ(original.rows(), 0);
    EXPECT_EQ(original.cols(), 0);
    EXPECT_EQ(original.size(), 0);
    EXPECT_EQ(original.data(), nullptr);

    // Verify moved array has correct values
    for (size_t i = 0; i < moved.rows(); ++i) {
        for (size_t j = 0; j < moved.cols(); ++j) {
            EXPECT_DOUBLE_EQ(moved(i, j), 5.0);
        }
    }
}

/**
 * @test AccessorBoundsCheck - Test boundary checking
 */
TEST(Array2DTest, AccessorBoundsCheck) {
    Array2D arr(3, 4);

    // Test valid access
    EXPECT_NO_THROW(arr(0, 0));
    EXPECT_NO_THROW(arr(2, 3));

    // Test row out of bounds
    EXPECT_THROW(arr(3, 0), std::out_of_range);
    EXPECT_THROW(arr(10, 0), std::out_of_range);

    // Test column out of bounds
    EXPECT_THROW(arr(0, 4), std::out_of_range);
    EXPECT_THROW(arr(0, 10), std::out_of_range);

    // Test both out of bounds
    EXPECT_THROW(arr(5, 5), std::out_of_range);

    // Test const version
    const Array2D& carr = arr;
    EXPECT_NO_THROW(carr(0, 0));
    EXPECT_THROW(carr(3, 0), std::out_of_range);
}

/**
 * @test Fill - Test fill method
 */
TEST(Array2DTest, Fill) {
    Array2D arr(4, 5);

    // Fill with positive value
    arr.fill(7.5);
    for (size_t i = 0; i < arr.rows(); ++i) {
        for (size_t j = 0; j < arr.cols(); ++j) {
            EXPECT_DOUBLE_EQ(arr(i, j), 7.5);
        }
    }

    // Fill with zero
    arr.fill(0.0);
    for (size_t i = 0; i < arr.rows(); ++i) {
        for (size_t j = 0; j < arr.cols(); ++j) {
            EXPECT_DOUBLE_EQ(arr(i, j), 0.0);
        }
    }

    // Fill with negative value
    arr.fill(-2.5);
    for (size_t i = 0; i < arr.rows(); ++i) {
        for (size_t j = 0; j < arr.cols(); ++j) {
            EXPECT_DOUBLE_EQ(arr(i, j), -2.5);
        }
    }
}

/**
 * @test CopyFrom - Test copy_from method
 */
TEST(Array2DTest, CopyFrom) {
    Array2D source(3, 4);
    double counter = 1.0;
    for (size_t i = 0; i < source.rows(); ++i) {
        for (size_t j = 0; j < source.cols(); ++j) {
            source(i, j) = counter++;
        }
    }

    Array2D dest(3, 4);

    // Test successful copy
    EXPECT_NO_THROW(dest.copy_from(source));
    for (size_t i = 0; i < dest.rows(); ++i) {
        for (size_t j = 0; j < dest.cols(); ++j) {
            EXPECT_DOUBLE_EQ(dest(i, j), source(i, j));
        }
    }

    // Test dimension mismatch
    Array2D wrong_size(4, 3);
    EXPECT_THROW(dest.copy_from(wrong_size), std::invalid_argument);

    Array2D wrong_rows(2, 4);
    EXPECT_THROW(dest.copy_from(wrong_rows), std::invalid_argument);

    Array2D wrong_cols(3, 5);
    EXPECT_THROW(dest.copy_from(wrong_cols), std::invalid_argument);
}

/**
 * @test Norms - Test various norm calculations
 */
TEST(Array2DTest, Norms) {
    Array2D arr(3, 3);

    // Fill with known values: 1, 2, 3, 4, 5, 6, 7, 8, 9
    double counter = 1.0;
    for (size_t i = 0; i < arr.rows(); ++i) {
        for (size_t j = 0; j < arr.cols(); ++j) {
            arr(i, j) = counter++;
        }
    }

    // Test max
    EXPECT_DOUBLE_EQ(arr.max(), 9.0);

    // Test min
    EXPECT_DOUBLE_EQ(arr.min(), 1.0);

    // Test L2 norm: sqrt(1^2 + 2^2 + ... + 9^2) = sqrt(285)
    double expected_l2 = std::sqrt(285.0);
    EXPECT_NEAR(arr.l2_norm(), expected_l2, 1e-10);

    // Test L-infinity norm: max(|element|) = 9
    EXPECT_DOUBLE_EQ(arr.linfty_norm(), 9.0);

    // Test with negative values
    Array2D arr_neg(2, 2);
    arr_neg(0, 0) = -5.0;
    arr_neg(0, 1) = 3.0;
    arr_neg(1, 0) = -2.0;
    arr_neg(1, 1) = 7.0;

    EXPECT_DOUBLE_EQ(arr_neg.max(), 7.0);
    EXPECT_DOUBLE_EQ(arr_neg.min(), -5.0);
    EXPECT_DOUBLE_EQ(arr_neg.linfty_norm(), 7.0);

    // L2 norm for negative array: sqrt(25 + 9 + 4 + 49) = sqrt(87)
    double expected_l2_neg = std::sqrt(87.0);
    EXPECT_NEAR(arr_neg.l2_norm(), expected_l2_neg, 1e-10);
}

/**
 * @test ArithmeticOperations - Test arithmetic operations
 */
TEST(Array2DTest, ArithmeticOperations) {
    // Create arrays for testing
    Array2D a(3, 3);
    Array2D b(3, 3);

    // Fill arrays
    double counter_a = 1.0;
    double counter_b = 10.0;
    for (size_t i = 0; i < a.rows(); ++i) {
        for (size_t j = 0; j < a.cols(); ++j) {
            a(i, j) = counter_a++;
            b(i, j) = counter_b++;
        }
    }

    // Test operator+=
    Array2D c = a;
    c += b;
    double counter_a_sum = 1.0;
    double counter_b_sum = 10.0;
    for (size_t i = 0; i < c.rows(); ++i) {
        for (size_t j = 0; j < c.cols(); ++j) {
            EXPECT_DOUBLE_EQ(c(i, j), counter_a_sum++ + counter_b_sum++);
        }
    }

    // Verify original arrays unchanged
    EXPECT_DOUBLE_EQ(a(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(b(0, 0), 10.0);

    // Test operator-=
    Array2D d = b;
    d -= a;
    double counter_a_diff = 1.0;
    double counter_b_diff = 10.0;
    for (size_t i = 0; i < d.rows(); ++i) {
        for (size_t j = 0; j < d.cols(); ++j) {
            EXPECT_DOUBLE_EQ(d(i, j), counter_b_diff++ - counter_a_diff++);
        }
    }

    // Test operator*=
    Array2D e = a;
    e *= 2.5;
    counter_a = 1.0;
    for (size_t i = 0; i < e.rows(); ++i) {
        for (size_t j = 0; j < e.cols(); ++j) {
            EXPECT_DOUBLE_EQ(e(i, j), counter_a++ * 2.5);
        }
    }

    // Test with negative scalar
    Array2D f = b;
    f *= -1.0;
    counter_b = 10.0;
    for (size_t i = 0; i < f.rows(); ++i) {
        for (size_t j = 0; j < f.cols(); ++j) {
            EXPECT_DOUBLE_EQ(f(i, j), -(counter_b++));
        }
    }

    // Test dimension mismatch for += and -=
    Array2D wrong_size(2, 3);
    EXPECT_THROW(a += wrong_size, std::invalid_argument);
    EXPECT_THROW(a -= wrong_size, std::invalid_argument);
}

/**
 * @test CopyAssignment - Test copy assignment operator
 */
TEST(Array2DTest, CopyAssignment) {
    Array2D a(3, 4, 5.0);
    Array2D b(2, 2, 1.0);

    b = a;  // This should resize b

    EXPECT_EQ(b.rows(), 3);
    EXPECT_EQ(b.cols(), 4);
    EXPECT_EQ(b.size(), 12);

    for (size_t i = 0; i < b.rows(); ++i) {
        for (size_t j = 0; j < b.cols(); ++j) {
            EXPECT_DOUBLE_EQ(b(i, j), 5.0);
        }
    }

    // Test self-assignment
    b = b;
    EXPECT_EQ(b.rows(), 3);
    EXPECT_EQ(b.cols(), 4);
}

/**
 * @test MoveAssignment - Test move assignment operator
 */
TEST(Array2DTest, MoveAssignment) {
    Array2D a(3, 4, 7.0);
    Array2D b(2, 2);

    b = std::move(a);

    EXPECT_EQ(b.rows(), 3);
    EXPECT_EQ(b.cols(), 4);
    EXPECT_EQ(b.size(), 12);

    for (size_t i = 0; i < b.rows(); ++i) {
        for (size_t j = 0; j < b.cols(); ++j) {
            EXPECT_DOUBLE_EQ(b(i, j), 7.0);
        }
    }

    EXPECT_EQ(a.rows(), 0);
    EXPECT_EQ(a.cols(), 0);
    EXPECT_EQ(a.size(), 0);

    // Test self-move (should be no-op)
    b = std::move(b);
    EXPECT_EQ(b.rows(), 3);
    EXPECT_EQ(b.cols(), 4);
}

/**
 * @test DataPointer - Test raw data pointer access
 */
TEST(Array2DTest, DataPointer) {
    Array2D arr(2, 3, 3.14);

    // Test non-const data()
    double* data = arr.data();
    EXPECT_NE(data, nullptr);
    EXPECT_DOUBLE_EQ(data[0], 3.14);

    // Modify through data pointer
    data[0] = 2.71;
    EXPECT_DOUBLE_EQ(arr(0, 0), 2.71);

    // Test const data()
    const Array2D& carr = arr;
    const double* cdata = carr.data();
    EXPECT_NE(cdata, nullptr);
    EXPECT_DOUBLE_EQ(cdata[0], 2.71);
}

/**
 * @test RowMajorLayout - Verify row-major memory layout
 */
TEST(Array2DTest, RowMajorLayout) {
    Array2D arr(2, 3);

    // Fill with sequential values
    double counter = 1.0;
    for (size_t i = 0; i < arr.rows(); ++i) {
        for (size_t j = 0; j < arr.cols(); ++j) {
            arr(i, j) = counter++;
        }
    }

    // Verify row-major layout in memory
    const double* data = arr.data();
    EXPECT_DOUBLE_EQ(data[0], 1.0);  // arr(0, 0)
    EXPECT_DOUBLE_EQ(data[1], 2.0);  // arr(0, 1)
    EXPECT_DOUBLE_EQ(data[2], 3.0);  // arr(0, 2)
    EXPECT_DOUBLE_EQ(data[3], 4.0);  // arr(1, 0)
    EXPECT_DOUBLE_EQ(data[4], 5.0);  // arr(1, 1)
    EXPECT_DOUBLE_EQ(data[5], 6.0);  // arr(1, 2)
}

/**
 * @test EdgeCases - Test edge cases and corner scenarios
 */
TEST(Array2DTest, EdgeCases) {
    // Test 1x1 array
    Array2D single(1, 1, 42.0);
    EXPECT_EQ(single.rows(), 1);
    EXPECT_EQ(single.cols(), 1);
    EXPECT_EQ(single.size(), 1);
    EXPECT_DOUBLE_EQ(single(0, 0), 42.0);
    EXPECT_DOUBLE_EQ(single.max(), 42.0);
    EXPECT_DOUBLE_EQ(single.min(), 42.0);
    EXPECT_DOUBLE_EQ(single.l2_norm(), 42.0);
    EXPECT_DOUBLE_EQ(single.linfty_norm(), 42.0);

    // Test 1xN array
    Array2D row_vec(1, 5);
    for (size_t j = 0; j < row_vec.cols(); ++j) {
        row_vec(0, j) = static_cast<double>(j);
    }
    EXPECT_EQ(row_vec.rows(), 1);
    EXPECT_EQ(row_vec.cols(), 5);
    EXPECT_DOUBLE_EQ(row_vec.max(), 4.0);
    EXPECT_DOUBLE_EQ(row_vec.min(), 0.0);

    // Test Nx1 array
    Array2D col_vec(5, 1);
    for (size_t i = 0; i < col_vec.rows(); ++i) {
        col_vec(i, 0) = static_cast<double>(i);
    }
    EXPECT_EQ(col_vec.rows(), 5);
    EXPECT_EQ(col_vec.cols(), 1);
    EXPECT_DOUBLE_EQ(col_vec.max(), 4.0);
    EXPECT_DOUBLE_EQ(col_vec.min(), 0.0);
}

/**
 * @test ChainedOperations - Test chaining of operations
 */
TEST(Array2DTest, ChainedOperations) {
    Array2D a(2, 2, 2.0);
    Array2D b(2, 2, 3.0);

    // Chain operations: (a + b) * 2
    Array2D c = a;
    (c += b) *= 2.0;

    for (size_t i = 0; i < c.rows(); ++i) {
        for (size_t j = 0; j < c.cols(); ++j) {
            EXPECT_DOUBLE_EQ(c(i, j), 10.0);  // (2 + 3) * 2 = 10
        }
    }
}

/**
 * @test LargeArray - Test with larger array sizes
 */
TEST(Array2DTest, LargeArray) {
    Array2D large(100, 100, 1.0);

    EXPECT_EQ(large.rows(), 100);
    EXPECT_EQ(large.cols(), 100);
    EXPECT_EQ(large.size(), 10000);

    // Verify corner elements
    EXPECT_DOUBLE_EQ(large(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(large(0, 99), 1.0);
    EXPECT_DOUBLE_EQ(large(99, 0), 1.0);
    EXPECT_DOUBLE_EQ(large(99, 99), 1.0);

    // Test norms
    EXPECT_DOUBLE_EQ(large.max(), 1.0);
    EXPECT_DOUBLE_EQ(large.min(), 1.0);
    EXPECT_DOUBLE_EQ(large.l2_norm(), 100.0);
    EXPECT_DOUBLE_EQ(large.linfty_norm(), 1.0);
}

// Main function for Google Test
