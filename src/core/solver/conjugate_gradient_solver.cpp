/**
 * @file conjugate_gradient_solver.cpp
 * @brief Implementation of Conjugate Gradient solver
 */

#include "conjugate_gradient_solver.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <chrono>
#include <algorithm>

/**
 * @brief Constructor for serial execution
 */
ConjugateGradientSolver::ConjugateGradientSolver(bool use_preconditioner, double lambda)
    : use_preconditioner_(use_preconditioner),
      lambda_(lambda),
      restart_threshold_(50),
      is_parallel_(false),
      comm_(MPI_COMM_NULL),
      rank_(0),
      size_(1),
      initialized_(false) {
    // Initialize neighbor ranks to -1 (no neighbors in serial mode)
    for (int i = 0; i < 4; ++i) {
        neighbor_rank_[i] = -1;
    }
}

/**
 * @brief Constructor for parallel execution
 */
ConjugateGradientSolver::ConjugateGradientSolver(bool use_preconditioner, double lambda,
                                                    MPI_Comm comm,
                                                    const int neighbor_rank[4],
                                                    int nx, int ny)
    : use_preconditioner_(use_preconditioner),
      lambda_(lambda),
      restart_threshold_(50),
      is_parallel_(true),
      comm_(comm),
      initialized_(false) {
    // Get MPI rank and size
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);

    // Copy neighbor ranks
    for (int i = 0; i < 4; ++i) {
        neighbor_rank_[i] = neighbor_rank[i];
    }
}

/**
 * @brief Destructor
 */
ConjugateGradientSolver::~ConjugateGradientSolver() = default;

/**
 * @brief Initialize work arrays
 */
void ConjugateGradientSolver::initialize_work_arrays(size_t rows, size_t cols) {
    r_ = std::make_unique<utils::Array2D>(rows, cols, 0.0);
    p_ = std::make_unique<utils::Array2D>(rows, cols, 0.0);
    Ap_ = std::make_unique<utils::Array2D>(rows, cols, 0.0);

    if (use_preconditioner_) {
        z_ = std::make_unique<utils::Array2D>(rows, cols, 0.0);
    }

    temp_ = std::make_unique<utils::Array2D>(rows, cols, 0.0);
    initialized_ = true;
}

/**
 * @brief Matrix-vector multiplication: y = Ax
 */
void ConjugateGradientSolver::matrix_vector_multiply(const utils::Array2D& x,
                                                     utils::Array2D& y) {
    // Exchange ghost cells first (required for parallel execution)
    if (is_parallel_) {
        exchange_ghost_cells(const_cast<utils::Array2D&>(x));
    }

    size_t rows = x.rows();
    size_t cols = x.cols();

    // Check dimensions
    if (y.rows() != rows || y.cols() != cols) {
        throw std::invalid_argument("Matrix-vector multiplication: dimension mismatch");
    }

    // Coefficient for diagonal
    double diag_coeff = 1.0 + 4.0 * lambda_;

    // For serial: compute on interior points
    // For parallel: compute on local interior (excluding ghost cells)
    size_t i_start = 1;
    size_t i_end = cols - 1;
    size_t j_start = 1;
    size_t j_end = rows - 1;

    // In serial mode, include boundaries
    if (!is_parallel_) {
        i_start = 0;
        i_end = cols;
        j_start = 0;
        j_end = rows;
    }

    // Compute y = Ax using 5-point stencil
    // y[i][j] = (1 + 4λ)x[i][j] - λ(x[i+1][j] + x[i-1][j] + x[i][j+1] + x[i][j-1])
    for (size_t j = j_start; j < j_end; ++j) {
        for (size_t i = i_start; i < i_end; ++i) {
            double north = (j > 0) ? x(j - 1, i) : 0.0;
            double south = (j < rows - 1) ? x(j + 1, i) : 0.0;
            double west = (i > 0) ? x(j, i - 1) : 0.0;
            double east = (i < cols - 1) ? x(j, i + 1) : 0.0;

            y(j, i) = diag_coeff *x(j, i) - lambda_ * (north + south + west + east);
        }
    }
}

/**
 * @brief Compute dot product: result = x^T y (with MPI reduction)
 */
double ConjugateGradientSolver::dot_product(const utils::Array2D& x,
                                              const utils::Array2D& y) {
    double local_dot = local_dot_product(x, y);

    if (is_parallel_) {
        double global_dot;
        MPI_Allreduce(&local_dot, &global_dot, 1, MPI_DOUBLE, MPI_SUM, comm_);
        return global_dot;
    }

    return local_dot;
}

/**
 * @brief Compute local dot product (no MPI reduction)
 */
double ConjugateGradientSolver::local_dot_product(const utils::Array2D& x,
                                                     const utils::Array2D& y) {
    double result = 0.0;
    size_t rows = x.rows();
    size_t cols = x.cols();

    for (size_t j = 0; j < rows; ++j) {
        for (size_t i = 0; i < cols; ++i) {
            result += x(j, i) * y(j, i);
        }
    }

    return result;
}

/**
 * @brief Apply Jacobi preconditioner: z = M^{-1} r
 */
void ConjugateGradientSolver::apply_preconditioner(const utils::Array2D& r,
                                                     utils::Array2D& z) {
    // For (I - λ∇²), the diagonal is (1 + 4λ)
    // So M^{-1} = 1/(1 + 4λ) * I
    double inv_diag = 1.0 / (1.0 + 4.0 * lambda_);

    size_t rows = r.rows();
    size_t cols = r.cols();

    for (size_t j = 0; j < rows; ++j) {
        for (size_t i = 0; i < cols; ++i) {
            z(j, i) = inv_diag * r(j, i);
        }
    }
}

/**
 * @brief Exchange ghost cells with neighbors
 */
void ConjugateGradientSolver::exchange_ghost_cells(utils::Array2D& x) {
    size_t rows = x.rows();
    size_t cols = x.cols();

    // Exchange with North neighbor (rank[0])
    if (neighbor_rank_[0] != -1) {
        MPI_Sendrecv(&x(0, 0), cols, MPI_DOUBLE, neighbor_rank_[0], 0,
                     &x(rows - 1, 0), cols, MPI_DOUBLE, neighbor_rank_[0], 0,
                     comm_, MPI_STATUS_IGNORE);
    }

    // Exchange with South neighbor (rank[1])
    if (neighbor_rank_[1] != -1) {
        MPI_Sendrecv(&x(rows - 2, 0), cols, MPI_DOUBLE, neighbor_rank_[1], 1,
                     &x(0, 0), cols, MPI_DOUBLE, neighbor_rank_[1], 1,
                     comm_, MPI_STATUS_IGNORE);
    }

    // Exchange with East neighbor (rank[2])
    if (neighbor_rank_[2] != -1) {
        // Need to send/receive column vectors
        std::vector<double> send_buffer(rows);
        std::vector<double> recv_buffer(rows);

        // Pack for send (last column before ghost)
        for (size_t j = 0; j < rows; j++) {
            send_buffer[j] = x(j, cols - 2);
        }

        MPI_Sendrecv(send_buffer.data(), rows, MPI_DOUBLE, neighbor_rank_[2], 2,
                     recv_buffer.data(), rows, MPI_DOUBLE, neighbor_rank_[2], 2,
                     comm_, MPI_STATUS_IGNORE);

        // Unpack to ghost column
        for (size_t j = 0; j < rows; j++) {
            x(j, cols - 1) = recv_buffer[j];
        }
    }

    // Exchange with West neighbor (rank[3])
    if (neighbor_rank_[3] != -1) {
        std::vector<double> send_buffer(rows);
        std::vector<double> recv_buffer(rows);

        // Pack for send (first column after ghost)
        for (size_t j = 0; j < rows; j++) {
            send_buffer[j] = x(j, 1);
        }

        MPI_Sendrecv(send_buffer.data(), rows, MPI_DOUBLE, neighbor_rank_[3], 3,
                     recv_buffer.data(), rows, MPI_DOUBLE, neighbor_rank_[3], 3,
                     comm_, MPI_STATUS_IGNORE);

        // Unpack to ghost column
        for (size_t j = 0; j < rows; j++) {
            x(j, 0) = recv_buffer[j];
        }
    }
}

/**
 * @brief Compute L2 norm (with MPI reduction)
 */
double ConjugateGradientSolver::norm(const utils::Array2D& x) {
    return std::sqrt(dot_product(x, x));
}

/**
 * @brief Compute local L2 norm (no MPI reduction)
 */
double ConjugateGradientSolver::local_norm(const utils::Array2D& x) {
    return std::sqrt(local_dot_product(x, x));
}

/**
 * @brief Print iteration progress
 */
void ConjugateGradientSolver::print_progress(int iteration, double residual) const {
    if (rank_ == 0 && iteration % 10 == 0) {
        std::cout << "  Iteration " << iteration
                  << ", Residual = " << std::scientific << residual << std::endl;
    }
}

/**
 * @brief Solve the linear system Ax = b
 */
void ConjugateGradientSolver::solve(const utils::Array2D& rhs,
                                     utils::Array2D& solution,
                                     const SolverParams& params) {
    // Update lambda from params if provided
    if (params.lambda > 0.0) {
        lambda_ = params.lambda;
    }

    // Check dimensions
    if (rhs.rows() != solution.rows() || rhs.cols() != solution.cols()) {
        throw std::invalid_argument("ConjugateGradientSolver: Dimension mismatch between rhs and solution");
    }

    size_t rows = rhs.rows();
    size_t cols = rhs.cols();

    // Initialize work arrays
    if (!initialized_ || r_->rows() != rows || r_->cols() != cols) {
        initialize_work_arrays(rows, cols);
    }

    // Reset statistics
    stats_.iterations = 0;
    stats_.restarts = 0;
    stats_.solve_time = 0.0;

    // Start timer
    auto start_time = std::chrono::high_resolution_clock::now();

    // Print solver info
    if (rank_ == 0 && params.verbose) {
        std::cout << "=== " << get_name() << " Solver ===" << std::endl;
        std::cout << "Grid size: " << rows << " x " << cols << std::endl;
        std::cout << "Lambda: " << lambda_ << std::endl;
        std::cout << "Tolerance: " << params.tolerance << std::endl;
        std::cout << "Max iterations: " << params.max_iterations << std::endl;
        std::cout << "Preconditioner: " << (use_preconditioner_ ? "Jacobi" : "None") << std::endl;
        std::cout << "Parallel: " << (is_parallel_ ? "Yes" : "No") << std::endl;
        if (is_parallel_) {
            std::cout << "Processes: " << size_ << std::endl;
        }
        std::cout << std::endl;
    }

    // Initialize: x = 0 (zero initial guess)
    solution.fill(0.0);

    // Compute initial residual: r = b - Ax = b (since x = 0)
    r_->copy_from(rhs);

    // Compute initial residual norm
    double r_norm = norm(*r_);

    // Check for trivial solution
    if (r_norm < params.tolerance) {
        stats_.residual = r_norm;
        stats_.convergence_rate = 0.0;

        auto end_time = std::chrono::high_resolution_clock::now();
        stats_.solve_time = std::chrono::duration<double>(end_time - start_time).count();

        if (rank_ == 0 && params.verbose) {
            std::cout << "Trivial solution (RHS norm below tolerance)" << std::endl;
            std::cout << "Final residual: " << std::scientific << stats_.residual << std::endl;
            std::cout << "Solve time: " << std::fixed << stats_.solve_time << " s" << std::endl;
        }

        return;
    }

    // Initial search direction
    double rTr = 0.0;  // Will be set correctly below
    if (use_preconditioner_) {
        // Preconditioned CG: z = M^{-1} r, p = z
        apply_preconditioner(*r_, *z_);
        p_->copy_from(*z_);
        rTr = dot_product(*r_, *z_);  // rTr = rᵀz = rᵀM^{-1}r
    } else {
        // Standard CG: p = r
        p_->copy_from(*r_);
        rTr = r_norm * r_norm;  // rTr = rᵀr
    }

    // Track previous residual for restart detection
    double previous_residual = r_norm;
    int iterations_without_progress = 0;

    // Main CG iteration loop
    for (size_t iteration = 1; iteration <= params.max_iterations; ++iteration) {
        stats_.iterations = iteration;

        // Compute Ap = A * p
        matrix_vector_multiply(*p_, *Ap_);

        // Compute alpha = (rᵀz) / (pᵀAp) for PCG, or (rᵀr) / (pᵀAp) for CG
        double pTAp = dot_product(*p_, *Ap_);

        // Check for breakdown
        if (std::abs(pTAp) < 1e-50) {
            if (rank_ == 0 && params.verbose) {
                std::cout << "Warning: Numerical breakdown at iteration " << iteration << std::endl;
                std::cout << "p^T A p = " << pTAp << std::endl;
            }
            break;
        }

        double alpha = rTr / pTAp;

        // Update solution: x = x + alpha * p
        for (size_t j = 0; j < rows; ++j) {
            for (size_t i = 0; i < cols; ++i) {
                solution(j, i) += alpha * (*p_)(j, i);
            }
        }

        // Compute new residual: r_new = r - alpha * Ap
        for (size_t j = 0; j < rows; ++j) {
            for (size_t i = 0; i < cols; ++i) {
                (*r_)(j, i) -= alpha * (*Ap_)(j, i);
            }
        }

        // Compute new residual norm
        r_norm = norm(*r_);

        // Check convergence
        if (r_norm < params.tolerance) {
            stats_.residual = r_norm;
            stats_.convergence_rate = std::pow(r_norm / (previous_residual + 1e-50),
                                               1.0 / static_cast<double>(iteration));

            auto end_time = std::chrono::high_resolution_clock::now();
            stats_.solve_time = std::chrono::duration<double>(end_time - start_time).count();

            if (rank_ == 0 && params.verbose) {
                std::cout << "Converged!" << std::endl;
                std::cout << "Iterations: " << stats_.iterations << std::endl;
                std::cout << "Restarts: " << stats_.restarts << std::endl;
                std::cout << "Final residual: " << std::scientific << stats_.residual << std::endl;
                std::cout << "Solve time: " << std::fixed << stats_.solve_time << " s" << std::endl;
                std::cout << "Avg time/iter: " << stats_.solve_time / stats_.iterations << " s" << std::endl;
            }

            print_progress(iteration, r_norm);
            return;
        }

        // Print progress
        if (params.verbose) {
            print_progress(iteration, r_norm);
        }

        // Check for restart (numerical stability)
        if (r_norm > previous_residual) {
            iterations_without_progress++;

            if (iterations_without_progress >= restart_threshold_) {
                stats_.restarts++;

                if (rank_ == 0 && params.verbose) {
                    std::cout << "Restarting at iteration " << iteration
                              << " (residual increased for " << iterations_without_progress
                              << " iterations)" << std::endl;
                }

                // Restart: recompute r = b - Ax with current x
                matrix_vector_multiply(solution, *temp_);
                for (size_t j = 0; j < rows; ++j) {
                    for (size_t i = 0; i < cols; ++i) {
                        (*r_)(j, i) = rhs(j, i) - (*temp_)(j, i);
                    }
                }
                r_norm = norm(*r_);

                // Reset search direction
                if (use_preconditioner_) {
                    apply_preconditioner(*r_, *z_);
                    p_->copy_from(*z_);
                    rTr = dot_product(*r_, *z_);  // rTr = rᵀz
                } else {
                    p_->copy_from(*r_);
                    rTr = r_norm * r_norm;  // rTr = rᵀr
                }

                iterations_without_progress = 0;
                previous_residual = r_norm;
                continue;
            }
        } else {
            iterations_without_progress = 0;
        }

        previous_residual = r_norm;

        // Compute beta = (r_newᵀz_new) / (rᵀz) for PCG, or (r_newᵀr_new) / (rᵀr) for CG
        double rTr_new = 0.0;
        if (use_preconditioner_) {
            // PCG: beta = (r_newᵀz_new) / (rᵀz)
            apply_preconditioner(*r_, *z_);
            rTr_new = dot_product(*r_, *z_);
        } else {
            // CG: beta = (r_newᵀr_new) / (rᵀr)
            rTr_new = r_norm * r_norm;
        }

        double beta = rTr_new / rTr;
        rTr = rTr_new;

        // Update search direction: p = z_new + beta * p (PCG) or p = r_new + beta * p (CG)
        if (use_preconditioner_) {
            // Preconditioned CG
            for (size_t j = 0; j < rows; ++j) {
                for (size_t i = 0; i < cols; ++i) {
                    (*p_)(j, i) = (*z_)(j, i) + beta * (*p_)(j, i);
                }
            }
        } else {
            // Standard CG
            for (size_t j = 0; j < rows; ++j) {
                for (size_t i = 0; i < cols; ++i) {
                    (*p_)(j, i) = (*r_)(j, i) + beta * (*p_)(j, i);
                }
            }
        }
    }

    // Did not converge within max iterations
    stats_.residual = r_norm;
    stats_.convergence_rate = std::pow(r_norm / (previous_residual + 1e-50),
                                       1.0 / static_cast<double>(stats_.iterations));

    auto end_time = std::chrono::high_resolution_clock::now();
    stats_.solve_time = std::chrono::duration<double>(end_time - start_time).count();

    if (rank_ == 0 && params.verbose) {
        std::cout << "Warning: Did not converge within " << params.max_iterations
                  << " iterations" << std::endl;
        std::cout << "Final residual: " << std::scientific << stats_.residual << std::endl;
        std::cout << "Solve time: " << std::fixed << stats_.solve_time << " s" << std::endl;
    }
}

/**
 * @brief Get solver name
 */
std::string ConjugateGradientSolver::get_name() const {
    if (use_preconditioner_) {
        return "Preconditioned Conjugate Gradient (PCG)";
    }
    return "Conjugate Gradient (CG)";
}

/**
 * @brief Reset solver state
 */
void ConjugateGradientSolver::reset() {
    stats_.iterations = 0;
    stats_.residual = 0.0;
    stats_.convergence_rate = 0.0;
    stats_.solve_time = 0.0;
    stats_.restarts = 0;
    initialized_ = false;
}
