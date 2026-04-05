/**
 * @file main.cpp
 * @brief Main program for 2D heat equation solver
 *
 * This program solves the 2D heat equation using implicit time integration
 * schemes and iterative linear solvers with MPI parallelization.
 *
 * Usage:
 *   heat_equation [options] [param_file]
 *
 * Example:
 *   heat_equation -f params.json
 *   heat_equation -nx 150 -ny 150 -nt 5 -stabp 0.1
 */

#include <iostream>
#include <memory>
#include <chrono>
#include <iomanip>

// MPI
#include "mpi/mpi_context.hpp"
#include "mpi/cartesian_topology.hpp"

// Mesh
#include "mesh/mesh2d.hpp"

// Solver
#include "core/solver/jacobi_solver.hpp"
#include "core/solver/sor_solver.hpp"
#include "core/solver/conjugate_gradient_solver.hpp"

// Time scheme
#include "core/scheme/implicit_euler.hpp"
#include "core/scheme/crank_nicolson.hpp"

// I/O
#include "io/param_reader.hpp"
#include "io/solution_exporter.hpp"

// Utilities
#include "utils/logger.hpp"

/**
 * @brief Create linear solver based on parameters
 * @param params Simulation parameters
 * @return Unique pointer to solver
 */
std::unique_ptr<ISolver> create_solver(const SimulationParams& params) {
    utils::Logger::get_instance().log_info(
        "Creating " + params.get_solver_type_string() + " solver"
    );

    SolverParams solver_params(
        params.solver_tolerance,
        params.max_solver_iterations,
        10,  // residual_check_interval
        true,  // compute_final_residual
        0.25  // lambda (relaxation parameter for SOR)
    );

    switch (params.solver_type) {
        case SolverType::Jacobi:
            return std::make_unique<JacobiSolver>();
        case SolverType::SOR:
            return std::make_unique<SORSolver>();
        case SolverType::ConjugateGradient:
            return std::make_unique<ConjugateGradientSolver>();
        default:
            utils::Logger::get_instance().log_warning(
                "Solver type not implemented, using Jacobi"
            );
            return std::make_unique<JacobiSolver>();
    }
}

/**
 * @brief Create time scheme based on parameters
 * @param params Simulation parameters
 * @return Unique pointer to time scheme
 */
std::unique_ptr<ITimeScheme> create_time_scheme(const SimulationParams& params) {
    utils::Logger::get_instance().log_info(
        "Creating " + params.get_scheme_type_string() + " time scheme"
    );

    switch (params.scheme_type) {
        case TimeSchemeType::ImplicitEuler:
            return std::make_unique<ImplicitEuler>();
        case TimeSchemeType::CrankNicolson:
            return std::make_unique<CrankNicolson>();
        default:
            utils::Logger::get_instance().log_warning(
                "Time scheme not implemented, using ImplicitEuler"
            );
            return std::make_unique<ImplicitEuler>();
    }
}

/**
 * @brief Print simulation statistics
 * @param total_time Total simulation time
 * @param params Simulation parameters
 */
void print_statistics(double total_time, const SimulationParams& params) {
    std::cout << "\n=== Simulation Statistics ===\n";
    std::cout << "Total time: " << std::fixed << std::setprecision(3)
              << total_time << " seconds\n";
    std::cout << "Average time per step: " << std::fixed << std::setprecision(6)
              << total_time / params.nt << " seconds\n";
    std::cout << "Grid points: " << params.nx * params.ny << "\n";
    std::cout << "Time steps: " << params.nt << "\n";
    std::cout << "=============================\n";
}

/**
 * @brief Main simulation loop
 * @param params Simulation parameters
 * @param mesh Current mesh
 * @param time_scheme Time integration scheme
 * @param exporter Solution exporter
 */
void run_simulation(const SimulationParams& params,
                    Mesh2D& mesh,
                    ITimeScheme& time_scheme,
                    SolutionExporter& exporter) {
    utils::Logger::get_instance().log_info("Starting simulation");

    // Time integration parameters
    TimeParams time_params(
        params.dt,
        mesh.hx(),
        mesh.hy(),
        params.alpha,
        params.t_start,
        params.t_end,
        params.nt
    );

    // Solver parameters
    SolverParams solver_params(
        params.solver_tolerance,
        params.max_solver_iterations
    );

    // Create next time step mesh
    Mesh2D next_mesh(mesh);

    // Export initial condition
    if (params.output_initial) {
        exporter.export_initial(mesh);
    }

    // Start timing
    auto start_time = std::chrono::high_resolution_clock::now();

    // Time integration loop
    for (size_t step = 0; step < params.nt; ++step) {
        utils::Logger::get_instance().log_info(
            "Time step " + std::to_string(step + 1) + "/" + std::to_string(params.nt)
        );

        // Perform time step
        time_scheme.step(mesh, next_mesh, time_params, solver_params);

        // Swap meshes
        std::swap(mesh, next_mesh);

        // Output if needed
        if ((step + 1) % params.output_interval == 0) {
            double current_time = params.t_start + (step + 1) * params.dt;
            exporter.export_solution(mesh, current_time, step + 1);
        }
    }

    // End timing
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    // Export final solution
    if (params.output_final) {
        exporter.export_final(mesh, params.t_end);
    }

    // Print statistics
    print_statistics(elapsed.count(), params);

    utils::Logger::get_instance().log_info("Simulation completed");
}

/**
 * @brief Main function
 * @param argc Argument count
 * @param argv Argument vector
 * @return Exit code
 */
int main(int argc, char** argv) {
    try {
        // Initialize MPI
        std::unique_ptr<MPIContext> mpi_context;
        std::unique_ptr<CartesianTopology> topology;

        // Check MPI support
        bool use_mpi = true;  // Default to MPI

        // Parse command line arguments first to check for -no-mpi flag
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "-no-mpi" || arg == "--no-mpi") {
                use_mpi = false;
                break;
            }
        }

        // Initialize MPI if enabled
        if (use_mpi) {
            mpi_context = std::make_unique<MPIContext>(argc, argv);
            topology = std::make_unique<CartesianTopology>();

            utils::Logger::get_instance().set_rank(mpi_context->rank());
            utils::Logger::get_instance().set_size(mpi_context->size());
            utils::Logger::get_instance().log_info(
                "MPI initialized: rank " + std::to_string(mpi_context->rank()) +
                " of " + std::to_string(mpi_context->size())
            );
        } else {
            utils::Logger::get_instance().log_info("Running in serial mode (no MPI)");
        }

        // Create parameter reader
        ParamReader param_reader;

        // Create default parameters
        SimulationParams params;

        // Parse command line
        std::string filename = param_reader.parse_command_line(argc, argv, params);

        // Load parameters from file if specified
        if (!filename.empty()) {
            params = param_reader.read_from_file(filename);

            // Re-apply command line overrides (command line takes precedence)
            // This is a simplified approach - a full implementation would
            // properly merge file and command line parameters
        }

        // Validate parameters
        params.validate();

        // Print parameters
        if (!use_mpi || mpi_context->is_root()) {
            ParamReader::print_params(params);
        }

        // Create mesh
        utils::Logger::get_instance().log_info("Creating mesh");
        Mesh2D mesh;

        if (use_mpi && topology) {
            mesh = Mesh2D(params.nx, params.ny, params.lx, params.ly, *topology);
        } else {
            mesh = Mesh2D(params.nx, params.ny, params.lx, params.ly);
        }

        // Initialize mesh (e.g., with initial condition)
        mesh.fill(0.0);  // Zero initial condition
        mesh.apply_dirichlet_bc(0.0, 0.0);  // Zero boundary conditions

        // Create time scheme
        auto time_scheme = create_time_scheme(params);

        // Create solution exporter
        bool mpi_mode = use_mpi && topology && mpi_context;
        SolutionExporter exporter(
            params.output_prefix,
            params.output_format,
            mpi_mode
        );

        // Run simulation
        run_simulation(params, mesh, *time_scheme, exporter);

        // MPI barrier to ensure all processes finish
        if (mpi_context) {
            mpi_context->barrier();
        }

        utils::Logger::get_instance().log_info("Program completed successfully");

        return 0;

    } catch (const std::exception& e) {
        utils::Logger::get_instance().log_error("Error: " + std::string(e.what()));
        return 1;
    } catch (...) {
        utils::Logger::get_instance().log_error("Unknown error occurred");
        return 1;
    }
}
