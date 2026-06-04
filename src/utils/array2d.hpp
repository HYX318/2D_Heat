/**
 * @file array2d.hpp
 * @brief RAII 2D array wrapper using modern C++17 features
 *
 * This header-only class provides a safe, efficient 2D array wrapper
 * with automatic memory management, bounds checking, and various utility methods.
 */

#ifndef ARRAY2D_HPP
#define ARRAY2D_HPP

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <algorithm>
#include <cmath>

namespace utils {

/**
 * @class Array2D
 * @brief RAII wrapper for 2D arrays with automatic memory management
 *
 * Features:
 * - Row-major contiguous memory layout
 * - Bounds-checked access via operator()
 * - Move semantics for efficient transfers
 * - Various utility methods (fill, copy, norms)
 * - Arithmetic operations
 * - Exception safety
 */
class Array2D {
public:
    /**
     * @brief Default constructor - creates uninitialized 2D array
     * @param rows Number of rows
     * @param cols Number of columns
     * @throws std::invalid_argument if rows or cols is 0
     */
    Array2D(size_t rows, size_t cols)
        : rows_(rows), cols_(cols), data_(std::make_unique<double[]>(rows * cols)) {
        if (rows == 0 || cols == 0) {
            throw std::invalid_argument("Array dimensions must be non-zero");
        }
    }

    /**
     * @brief Fill constructor - creates 2D array filled with specified value
     * @param rows Number of rows
     * @param cols Number of columns
     * @param value Fill value for all elements
     * @throws std::invalid_argument if rows or cols is 0
     */
    Array2D(size_t rows, size_t cols, double value)
        : Array2D(rows, cols) {
        fill(value);
    }

    /**
     * @brief Copy constructor - performs deep copy
     * @param other Array to copy from
     */
    Array2D(const Array2D& other)
        : rows_(other.rows_), cols_(other.cols_),
          data_(std::make_unique<double[]>(other.rows_ * other.cols_)) {
        std::copy(other.data_.get(), other.data_.get() + other.rows_ * other.cols_,
                  data_.get());
    }

    /**
     * @brief Move constructor - efficient transfer of ownership
     * @param other Array to move from
     */
    Array2D(Array2D&& other) noexcept
        : rows_(other.rows_), cols_(other.cols_), data_(std::move(other.data_)) {
        other.rows_ = 0;
        other.cols_ = 0;
    }

    /**
     * @brief Destructor - automatically handled by unique_ptr
     */
    ~Array2D() = default;

    /**
     * @brief Copy assignment operator
     * @param other Array to copy from
     * @return Reference to this array
     */
    Array2D& operator=(const Array2D& other) {
        if (this != &other) {
            Array2D temp(other);
            swap(temp);
        }
        return *this;
    }

    /**
     * @brief Move assignment operator
     * @param other Array to move from
     * @return Reference to this array
     */
    Array2D& operator=(Array2D&& other) noexcept {
        if (this != &other) {
            rows_ = other.rows_;
            cols_ = other.cols_;
            data_ = std::move(other.data_);
            other.rows_ = 0;
            other.cols_ = 0;
        }
        return *this;
    }

    /**
     * @brief Bounds-checked element access (mutable)
     * @param i Row index (0-based)
     * @param j Column index (0-based)
     * @return Reference to element at (i, j)
     * @throws std::out_of_range if indices are out of bounds
     */
    double& operator()(size_t i, size_t j) {
        check_bounds(i, j);
        return data_[i * cols_ + j];
    }

    /**
     * @brief Bounds-checked element access (const)
     * @param i Row index (0-based)
     * @param j Column index (0-based)
     * @return Const reference to element at (i, j)
     * @throws std::out_of_range if indices are out of bounds
     */
    const double& operator()(size_t i, size_t j) const {
        check_bounds(i, j);
        return data_[i * cols_ + j];
    }

    /**
     * @brief Fast element access without bounds checking.
     *
     * Intended for internal numerical kernels that already validate loop bounds.
     */
    double& unchecked(size_t i, size_t j) noexcept {
        return data_[i * cols_ + j];
    }

    /**
     * @brief Fast const element access without bounds checking.
     */
    const double& unchecked(size_t i, size_t j) const noexcept {
        return data_[i * cols_ + j];
    }

    /**
     * @brief Pointer to the first element of a row.
     */
    double* row_data(size_t i) noexcept {
        return data_.get() + i * cols_;
    }

    /**
     * @brief Const pointer to the first element of a row.
     */
    const double* row_data(size_t i) const noexcept {
        return data_.get() + i * cols_;
    }

    /**
     * @brief Get number of rows
     * @return Number of rows
     */
    size_t rows() const noexcept {
        return rows_;
    }

    /**
     * @brief Get number of columns
     * @return Number of columns
     */
    size_t cols() const noexcept {
        return cols_;
    }

    /**
     * @brief Get total number of elements
     * @return rows * cols
     */
    size_t size() const noexcept {
        return rows_ * cols_;
    }

    /**
     * @brief Fill all elements with specified value
     * @param value Value to fill with
     */
    void fill(double value) noexcept {
        std::fill(data_.get(), data_.get() + rows_ * cols_, value);
    }

    /**
     * @brief Copy data from another array
     * @param other Source array
     * @throws std::invalid_argument if dimensions don't match
     */
    void copy_from(const Array2D& other) {
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            throw std::invalid_argument("Array dimensions must match for copy_from");
        }
        std::copy(other.data_.get(), other.data_.get() + rows_ * cols_,
                  data_.get());
    }

    /**
     * @brief Compute maximum element value
     * @return Maximum value in the array
     */
    double max() const noexcept {
        return *std::max_element(data_.get(), data_.get() + rows_ * cols_);
    }

    /**
     * @brief Compute minimum element value
     * @return Minimum value in the array
    */
    double min() const noexcept {
        return *std::min_element(data_.get(), data_.get() + rows_ * cols_);
    }

    /**
     * @brief Compute L2 norm (Euclidean norm)
     * @return sqrt(sum of squares of all elements)
     */
    double l2_norm() const noexcept {
        double sum = 0.0;
        for (size_t i = 0; i < rows_ * cols_; ++i) {
            sum += data_[i] * data_[i];
        }
        return std::sqrt(sum);
    }

    /**
     * @brief Compute L-infinity norm (maximum absolute value)
     * @return Maximum absolute value in the array
     */
    double linfty_norm() const noexcept {
        double max_abs = 0.0;
        for (size_t i = 0; i < rows_ * cols_; ++i) {
            max_abs = std::max(max_abs, std::abs(data_[i]));
        }
        return max_abs;
    }

    /**
     * @brief In-place addition with another array
     * @param other Array to add
     * @return Reference to this array
     * @throws std::invalid_argument if dimensions don't match
     */
    Array2D& operator+=(const Array2D& other) {
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            throw std::invalid_argument("Array dimensions must match for addition");
        }
        for (size_t i = 0; i < rows_ * cols_; ++i) {
            data_[i] += other.data_[i];
        }
        return *this;
    }

    /**
     * @brief In-place subtraction with another array
     * @param other Array to subtract
     * @return Reference to this array
     * @throws std::invalid_argument if dimensions don't match
     */
    Array2D& operator-=(const Array2D& other) {
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            throw std::invalid_argument("Array dimensions must match for subtraction");
        }
        for (size_t i = 0; i < rows_ * cols_; ++i) {
            data_[i] -= other.data_[i];
        }
        return *this;
    }

    /**
     * @brief In-place scalar multiplication
     * @param scalar Scalar to multiply with
     * @return Reference to this array
     */
    Array2D& operator*=(double scalar) noexcept {
        for (size_t i = 0; i < rows_ * cols_; ++i) {
            data_[i] *= scalar;
        }
        return *this;
    }

    /**
     * @brief Get raw pointer to underlying data
     * @return Pointer to contiguous memory
     */
    double* data() noexcept {
        return data_.get();
    }

    /**
     * @brief Get raw pointer to underlying data (const)
     * @return Const pointer to contiguous memory
     */
    const double* data() const noexcept {
        return data_.get();
    }

private:
    /**
     * @brief Check if indices are within bounds
     * @param i Row index
     * @param j Column index
     * @throws std::out_of_range if indices are invalid
     */
    void check_bounds(size_t i, size_t j) const {
        if (i >= rows_) {
            throw std::out_of_range("Row index out of range");
        }
        if (j >= cols_) {
            throw std::out_of_range("Column index out of range");
        }
    }

    /**
     * @brief Swap contents with another array
     * @param other Array to swap with
     */
    void swap(Array2D& other) noexcept {
        std::swap(rows_, other.rows_);
        std::swap(cols_, other.cols_);
        std::swap(data_, other.data_);
    }

    size_t rows_;                          ///< Number of rows
    size_t cols_;                          ///< Number of columns
    std::unique_ptr<double[]> data_;       ///< Contiguous data storage
};

} // namespace utils

#endif // ARRAY2D_HPP
