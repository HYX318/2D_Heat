/**
 * @file sor_solver.cpp
 * @brief Implementation of SOR solver
 */

#include "sor_solver.hpp"
#include "../../mpi/ghost_cell_exchange.hpp"
#include <algorithm>
#include <iostream>

// Serial constructor
SORSolver::SORSolver(double omega)
    : omega_(omega)
    , auto_omega_(1.0)
    , use_auto_omega_(omega == 0.0)
    , use_red_black_(false)
    , is_parallel_(false)
    , ghost_exchange_(nullptr)
    , comm_(MPI_COMM_NULL)
    , mpi_rank_(0)
    , mpi_size_(1)
{
    validate_omega(omega);
}

// Parallel constructor
SORSolver::SORSolver(double omega, GhostCellExchange* ghost_exchange, MPI_Comm comm)
    : omega_(omega)
    , auto_omega_(1.0)
    , use_auto_omega_(omega == 0.0)
    , use_red_black_(true)  // Use red-black by default in parallel
    , is_parallel_(true)
    , ghost_exchange_(ghost_exchange)
    , comm_(comm)
    , mpi_rank_(0)
    , mpi_size_(1)
{
    validate_omega(omega);
    init_mpi();

    if (ghost_exchange_ == nullptr) {
        throw std::invalid_argument("GhostCellExchange cannot be null in parallel mode");
    }
}

// Move constructor
SORSolver::SORSolver(SORSolver&& other) noexcept
    : omega_(other.omega_)
    , auto_omega_(other.auto_omega_)
    , use_auto_omega_(other.use_auto_omega_)
    , use_red_black_(other.use_red_black_)
    , is_parallel_(other.is_parallel_)
    , ghost_exchange_(other.ghost_exchange_)
    , comm_(other.comm_)
    , mpi_rank_(other.mpi_rank_)
    , mpi_size_(other.mpi_size_)
    , stats_(other.stats_)
{
    other.ghost_exchange_ = nullptr;
    other.comm_ = MPI_COMM_NULL;
    other.is_parallel_ = false;
}

// Move assignment
SORSolver& SORSolver::operator=(SORSolver&& other) noexcept {
    if (this != &other) {
        omega_ = other.omega_;
        auto_omega_ = other.auto_omega_;
        use_auto_omega_ = other.use_auto_omega_;
        use_red_black_ = other.use_red_black_;
        is_parallel_ = other.is_parallel_;
        ghost_exchange_ = other.ghost_exchange_;
        comm_ = other.comm_;
        mpi_rank_ = other.mpi_rank_;
        mpi_size_ = other.mpi_size_;
        stats_ = other.stats_;

        other.ghost_exchange_ = nullptr;
        other.comm_ = MPI_COMM_NULL;
        other.is_parallel_ = false;
    }
    return *this;
}

// Main solve method
void SORSolver::solve(const utils::Array2D& rhs, utils::Array2D& solution,
                       const SolverParams& params) {
    // Reset statistics
    stats_.reset();
    auto start_time = std::chrono::high_resolution_clock::now();

    // Validate arrays
    validate_arrays(rhs, solution);

    // Check if we're actually in parallel mode
    bool actually_parallel = is_parallel_ && check_parallel();

    // Compute optimal omega if needed
    if (use_auto_omega_) {
        auto_omega_ = compute_optimal_omega(rhs.cols(), rhs.rows());
        omega_ = auto_omega_;
    }

    // Compute initial residual
    double residual = compute_residual(rhs, solution, params.lambda);
    stats_.initial_residual = residual;

    // Check if already converged
    if (residual < params.tolerance) {
        stats_.converged = true;
        stats_.iterations = 0;
        stats_.final_residual = residual;
        stats_.reduction_factor = 1.0;
        auto end_time = std::chrono::high_resolution_clock::now();
        stats_.solve_time = std::chrono::duration<double>(end_time - start_time).count();
        return;
    }

    // SOR iteration loop
    size_t iter;
    for (iter = 0; iter < params.max_iterations; ++iter) {
        // Perform one SOR iteration
        if (actually_parallel) {
            sor_iteration_parallel(rhs, solution, params.lambda);
        } else if (use_red_black_) {
            // Red-black ordering: update red points, then black points
            sor_iteration_redblack(rhs, solution, params.lambda, 0);  // Red
            sor_iteration_redblack(rhs, solution, params.lambda, 1);  // Black
        } else {
            sor_iteration_dict(rhs, solution, params.lambda);
        }

        // Compute residual at specified intervals
        if (iter % params.residual_check_interval == 0 || iter == params.max_iterations - 1) {
            residual = compute_residual(rhs, solution, params.lambda);
            
            // Check convergence
            if (residual < params.tolerance) {
                stats_.converged = true;
                break;
            }
        }
    }

    // Final statistics
    stats_.iterations = iter + 1;
    stats_.final_residual = residual;
    stats_.reduction_factor = residual / stats_.initial_residual;

    // Compute final residual if requested
    if (params.compute_final_residual) {
        stats_.final_residual = compute_residual(rhs, solution, params.lambda);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    stats_.solve_time = std::chrono::duration<double>(end_time - start_time).count();
}

// Get statistics
SolverStats SORSolver::get_stats() const {
    return stats_;
}

// Get solver name
std::string SORSolver::get_name() const {
    std::string name = "SOR";
    if (use_red_black_) {
        name += "-RedBlack";
    }
    if (is_parallel_ && check_parallel()) {
        name += "-Parallel";
    }
    return name;
}

// Get solver type
SolverType SORSolver::get_type() const {
    return SolverType::SOR;
}

// Reset solver state
void SORSolver::reset() {
    stats_.reset();
}

// Set omega
void SORSolver::set_omega(double omega) {
    validate_omega(omega);
    omega_ = omega;
    use_auto_omega_ = (omega == 0.0);
}

// Get omega
double SORSolver::get_omega() const {
    return omega_;
}

// Check if auto omega
bool SORSolver::is_auto_omega() const {
    return use_auto_omega_;
}

// Enable red-black ordering
void SORSolver::enable_red_black(bool enable) {
    use_red_black_ = enable;
}

// Check red-black ordering
bool SORSolver::is_red_black() const {
    return use_red_black_;
}

// Compute optimal omega
double SORSolver::compute_optimal_omega(size_t nx, size_t ny) {
    // For the discrete Poisson equation, spectral radius of Jacobi iteration matrix:
    // rho = cos(pi/(nx+1)) + cos(pi/(ny+1))
    double rho = std::cos(M_PI / (nx + 1)) + std::cos(M_PI / (ny + 1));
    
    // Optimal relaxation parameter:
    // omega_opt = 2 / (1 + sqrt(1 - rho^2))
    double omega_opt = 2.0 / (1.0 + std::sqrt(1.0 - rho * rho));
    
    return omega_opt;
}

// Validate omega
void SORSolver::validate_omega(double omega) {
    if (omega >= 2.0 && omega != 0.0) {
        throw std::invalid_argument("Omega must be in range [0, 2) for convergence");
    }
    if (omega < 0.0) {
        throw std::invalid_argument("Omega must be non-negative");
    }
}

// Compute residual
double SORSolver::compute_residual(const utils::Array2D& rhs,
                                    const utils::Array2D& solution,
                                    double lambda) {
    double max_residual = 0.0;
    size_t nx = solution.cols();
    size_t ny = solution.rows();

    // Compute residual at interior points (excluding boundaries)
    for (size_t i = 1; i < ny - 1; ++i) {
        for (size_t j = 1; j < nx - 1; ++j) {
            // Residual = Ax - b (for Poisson: -4*x + neighbors - b)
            // For the discretized heat equation:
            // residual = rhs[i][j] - solution[i][j] +
            //             lambda * (solution[i+1][j] + solution[i-1][j] +
            //                      solution[i][j+1] + solution[i][j-1] - 4*solution[i][j])
            
            double laplacian = solution(i+1, j) + solution(i-1, j) +
                               solution(i, j+1) + solution(i, j-1) - 4.0 * solution(i, j);
            
            double residual = std::abs(rhs(i, j) - solution(i, j) + lambda * laplacian);
            max_residual = std::max(max_residual, residual);
        }
    }

    return max_residual;
}

// SOR iteration with dictionary ordering
void SORSolver::sor_iteration_dict(const utils::Array2D& rhs,
                                    utils::Array2D& solution,
                                    double lambda) {
    size_t nx = solution.cols();
    size_t ny = solution.rows();
    double coeff = 1.0 / (1.0 - 4.0 * lambda);

    // Update in lexicographic order
    for (size_t i = 1; i < ny - 1; ++i) {
        for (size_t j = 1; j < nx - 1; ++j) {
            double x_old = solution(i, j);
            
            // SOR update formula
            double sum = solution(i+1, j) + solution(i-1, j) +
                         solution(i, j+1) + solution(i, j-1);
            
            double x_new = (1.0 - omega_) * x_old + 
                           omega_ * coeff * (rhs(i, j) + lambda * sum);
            
            solution(i, j) = x_new;
        }
    }
}

// SOR iteration with red-black ordering
void SORSolver::sor_iteration_redblack(const utils::Array2D& rhs,
                                        utils::Array2D& solution,
                                        double lambda,
                                        int color) {
    size_t nx = solution.cols();
    size_t ny = solution.rows();
    double coeff = 1.0 / (1.0 - 4.0 * lambda);

    // Update points of specified color
    for (size_t i = 1; i < ny - 1; ++i) {
        for (size_t j = 1; j < nx - 1; ++j) {
            // Check if this point matches the color
            if ((i + j) % 2 != color) {
                continue;
            }

            double x_old = solution(i, j);
            
            double sum = solution(i+1, j) + solution(i-1, j) +
                         solution(i, j+1) + solution(i, j-1);
            
            double x_new = (1.0 - omega_) * x_old + 
                           omega_ * coeff * (rhs(i, j) + lambda * sum);
            
            solution(i, j) = x_new;
        }
    }
}

// Parallel SOR iteration with ghost cell exchange
void SORSolver::sor_iteration_parallel(const utils::Array2D& rhs,
                                        utils::Array2D& solution,
                                        double lambda) {
    if (ghost_exchange_ == nullptr) {
        throw std::runtime_error("GhostCellExchange is null in parallel mode");
    }

    // Exchange ghost cells before iteration
    ghost_exchange_->exchange(solution);

    size_t nx = solution.cols();
    size_t ny = solution.rows();
    double coeff = 1.0 / (1.0 - 4.0 * lambda);

    if (use_red_black_) {
        // Red-black ordering with ghost cells
        // Note: Interior points are at indices [1..ny] x [1..nx]
        // Ghost cells are at index 0 and ny+1, nx+1
        
        // Update red points
        for (size_t i = 1; i <= ny; ++i) {
            for (size_t j = 1; j <= nx; ++j) {
                if ((i + j) % 2 != 0) continue;  // Only red points
                
                double x_old = solution(i, j);
                double sum = solution(i+1, j) + solution(i-1, j) +
                             solution(i, j+1) + solution(i, j-1);
                double x_new = (1.0 - omega_) * x_old +
                               omega_ * coeff * (rhs(i, j) + lambda * sum);
                solution(i, j) = x_new;
            }
        }

        // Exchange ghost cells again
        ghost_exchange_->exchange(solution);

        // Update black points
        for (size_t i = 1; i <= ny; ++i) {
            for (size_t j = 1; j <= nx; ++j) {
                if ((i + j) % 2 != 1) continue;  // Only black points
                
                double x_old = solution(i, j);
                double sum = solution(i+1, j) + solution(i-1, j) +
                             solution(i, j+1) + solution(i, j-1);
                double x_new = (1.0 - omega_) * x_old +
                               omega_ * coeff * (rhs(i, j) + lambda * sum);
                solution(i, j) = x_new;
            }
        }
    } else {
        // Dictionary ordering (not recommended for parallel)
        for (size_t i = 1; i <= ny; ++i) {
            for (size_t j = 1; j <= nx; ++j) {
                double x_old = solution(i, j);
                double sum = solution(i+1, j) + solution(i-1, j) +
                             solution(i, j+1) + solution(i, j-1);
                double x_new = (1.0 - omega_) * x_old +
                               omega_ * coeff * (rhs(i, j) + lambda * sum);
                solution(i, j) = x_new;
            }
        }
    }
}

// Validate arrays
void SORSolver::validate_arrays(const utils::Array2D& rhs,
                                 const utils::Array2D& solution) {
    if (rhs.rows() != solution.rows() || rhs.cols() != solution.cols()) {
        throw std::invalid_argument("RHS and solution arrays must have the same dimensions");
    }
    
    if (rhs.rows() < 3 || rhs.cols() < 3) {
        throw std::invalid_argument("Arrays must be at least 3x3 (for interior points)");
    }
}

// Initialize MPI state
void SORSolver::init_mpi() {
    int initialized = 0;
    MPI_Initialized(&initialized);
    
    if (initialized) {
        MPI_Comm_rank(comm_, &mpi_rank_);
        MPI_Comm_size(comm_, &mpi_size_);
    } else {
        mpi_rank_ = 0;
        mpi_size_ = 1;
        is_parallel_ = false;
    }
}

// Check if running in parallel
bool SORSolver::check_parallel() const {
    int initialized = 0;
    MPI_Initialized(&initialized);
    return initialized && mpi_size_ > 1;
}
