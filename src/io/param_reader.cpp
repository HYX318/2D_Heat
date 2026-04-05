/**
 * @file param_reader.cpp
 * @brief Implementation of parameter reader
 */

#include "param_reader.hpp"
#include "../utils/logger.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>

// Validate parameters
void SimulationParams::validate() const {
    if (nx == 0 || ny == 0) {
        throw std::invalid_argument("Grid dimensions must be positive");
    }
    if (nx < 2 || ny < 2) {
        throw std::invalid_argument("Grid dimensions must be at least 2");
    }
    if (lx <= 0.0 || ly <= 0.0) {
        throw std::invalid_argument("Physical dimensions must be positive");
    }
    if (dt <= 0.0) {
        throw std::invalid_argument("Time step must be positive");
    }
    if (t_start >= t_end) {
        throw std::invalid_argument("End time must be greater than start time");
    }
    if (nt == 0) {
        throw std::invalid_argument("Number of time steps must be positive");
    }
    if (alpha <= 0.0) {
        throw std::invalid_argument("Diffusion coefficient must be positive");
    }
    if (solver_tolerance <= 0.0) {
        throw std::invalid_argument("Solver tolerance must be positive");
    }
    if (max_solver_iterations == 0) {
        throw std::invalid_argument("Maximum solver iterations must be positive");
    }
    if (output_interval == 0) {
        throw std::invalid_argument("Output interval must be positive");
    }
}

double SimulationParams::compute_dt() const {
    // Compute time step from stability parameter
    // For explicit schemes: dt = stabp * dx^2 / alpha
    // For implicit schemes, this is just a guideline
    double dx = lx / (nx - 1);
    double dy = ly / (ny - 1);
    double min_d2 = std::min(dx * dx, dy * dy);
    return stabp * min_d2 / alpha;
}

double SimulationParams::compute_stabp() const {
    // Compute stability parameter from time step
    double dx = lx / (nx - 1);
    double dy = ly / (ny - 1);
    double min_d2 = std::min(dx * dx, dy * dy);
    return alpha * dt / min_d2;
}

std::string SimulationParams::get_solver_type_string() const {
    switch (solver_type) {
        case SolverType::Jacobi: return "Jacobi";
        case SolverType::GaussSeidel: return "GaussSeidel";
        case SolverType::SOR: return "SOR";
        case SolverType::ConjugateGradient: return "ConjugateGradient";
        case SolverType::Multigrid: return "Multigrid";
        default: return "Unknown";
    }
}

std::string SimulationParams::get_scheme_type_string() const {
    switch (scheme_type) {
        case TimeSchemeType::ImplicitEuler: return "ImplicitEuler";
        case TimeSchemeType::CrankNicolson: return "CrankNicolson";
        case TimeSchemeType::ExplicitEuler: return "ExplicitEuler";
        case TimeSchemeType::RungeKutta4: return "RungeKutta4";
        default: return "Unknown";
    }
}

std::string SimulationParams::get_output_format_string() const {
    switch (output_format) {
        case OutputFormat::Text: return "Text";
        case OutputFormat::VTK: return "VTK";
        case OutputFormat::HDF5: return "HDF5";
        default: return "Unknown";
    }
}

ParamReader::ParamReader() {}

InputFormat ParamReader::detect_format(const std::string& filename) const {
    // Check file extension
    if (filename.size() > 5) {
        std::string ext = filename.substr(filename.size() - 5);
        if (ext == ".json" || ext == ".JSON") {
            return InputFormat::JSON;
        }
    }
    // Default to text format
    return InputFormat::Text;
}

std::string ParamReader::trim(const std::string& str) const {
    size_t start = str.find_first_not_of(" \t\n\r");
    if (start == std::string::npos) return "";
    size_t end = str.find_last_not_of(" \t\n\r");
    return str.substr(start, end - start + 1);
}

std::vector<std::string> ParamReader::split(const std::string& str,
                                              char delim) const {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delim)) {
        std::string trimmed = trim(token);
        if (!trimmed.empty()) {
            tokens.push_back(trimmed);
        }
    }
    return tokens;
}

SolverType ParamReader::parse_solver_type(const std::string& str) const {
    std::string lower = str;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c){ return std::tolower(c); });

    if (lower == "jacobi" || lower == "jacobi_solver") {
        return SolverType::Jacobi;
    } else if (lower == "gaussseidel" || lower == "gauss_seidel") {
        return SolverType::GaussSeidel;
    } else if (lower == "sor" || lower == "sor_solver") {
        return SolverType::SOR;
    } else if (lower == "conjugategradient" || lower == "cg" ||
               lower == "conjugate_gradient") {
        return SolverType::ConjugateGradient;
    } else if (lower == "multigrid" || lower == "mg") {
        return SolverType::Multigrid;
    } else {
        throw std::invalid_argument("Unknown solver type: " + str);
    }
}

TimeSchemeType ParamReader::parse_scheme_type(const std::string& str) const {
    std::string lower = str;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c){ return std::tolower(c); });

    if (lower == "impliciteuler" || lower == "implicit_euler" ||
        lower == "backward_euler" || lower == "be") {
        return TimeSchemeType::ImplicitEuler;
    } else if (lower == "cranknicolson" || lower == "crank_nicolson" ||
               lower == "cn") {
        return TimeSchemeType::CrankNicolson;
    } else if (lower == "expliciteuler" || lower == "explicit_euler" ||
               lower == "fe") {
        return TimeSchemeType::ExplicitEuler;
    } else if (lower == "rungekutta4" || lower == "runge_kutta4" ||
               lower == "rk4") {
        return TimeSchemeType::RungeKutta4;
    } else {
        throw std::invalid_argument("Unknown time scheme type: " + str);
    }
}

OutputFormat ParamReader::parse_output_format(const std::string& str) const {
    std::string lower = str;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c){ return std::tolower(c); });

    if (lower == "text" || lower == "txt") {
        return OutputFormat::Text;
    } else if (lower == "vtk") {
        return OutputFormat::VTK;
    } else if (lower == "hdf5" || lower == "h5") {
        return OutputFormat::HDF5;
    } else {
        throw std::invalid_argument("Unknown output format: " + str);
    }
}

SimulationParams ParamReader::read_from_file(const std::string& filename,
                                               InputFormat format) {
    // Auto-detect format if needed
    if (format == InputFormat::Auto) {
        format = detect_format(filename);
    }

    if (format == InputFormat::JSON) {
        return read_json_format(filename);
    } else {
        return read_text_format(filename);
    }
}

SimulationParams ParamReader::read_text_format(const std::string& filename) {
    utils::Logger::get_instance().log_info(
        "Reading parameters from text file: " + filename
    );

    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    SimulationParams params;
    std::string line;
    std::vector<std::string> values;

    // Read all values from file (backward compatible format)
    while (std::getline(file, line)) {
        // Skip empty lines and comments
        std::string trimmed = trim(line);
        if (trimmed.empty() || trimmed[0] == '#') {
            continue;
        }
        values.push_back(trimmed);
    }

    if (values.empty()) {
        throw std::runtime_error("No data found in file: " + filename);
    }

    try {
        // Parse values (backward compatible with old Param.in format)
        // Format: nx ny (line 1)
        //          nt (line 2)
        //          stabp (line 3)

        if (values.size() >= 1) {
            auto tokens = split(values[0], ' ');
            if (tokens.size() >= 2) {
                params.nx = std::stoul(tokens[0]);
                params.ny = std::stoul(tokens[1]);
            }
        }

        if (values.size() >= 2) {
            params.nt = std::stoul(values[1]);
        }

        if (values.size() >= 3) {
            params.stabp = std::stod(values[2]);
        }

        // Compute derived parameters
        params.dt = params.compute_dt();
        params.t_end = params.nt * params.dt;

        // Validate parameters
        params.validate();

    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing text file: " + std::string(e.what()));
    }

    utils::Logger::get_instance().log_info("Parameters loaded successfully");
    return params;
}

SimulationParams ParamReader::read_json_format(const std::string& filename) {
    utils::Logger::get_instance().log_info(
        "Reading parameters from JSON file: " + filename
    );

    // Note: This is a simplified JSON parser
    // For production code, consider using a proper JSON library like nlohmann/json

    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    SimulationParams params;
    std::string content((std::istreambuf_iterator<char>(file)),
                         std::istreambuf_iterator<char>());

    try {
        // Very simple JSON parsing (not production quality)
        // In practice, use nlohmann/json or similar library

        // For now, just read as text format for backward compatibility
        utils::Logger::get_instance().log_warning(
            "JSON parsing not fully implemented, falling back to text format"
        );
        return read_text_format(filename);

    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing JSON file: " + std::string(e.what()));
    }
}

std::string ParamReader::parse_command_line(int argc, char** argv,
                                              SimulationParams& params) {
    std::string filename;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            print_help(argv[0]);
            throw std::runtime_error("Help requested");
        }
        else if (arg == "-f" || arg == "--file") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing filename after " + arg);
            }
            filename = argv[++i];
        }
        else if (arg == "-nx" || arg == "--nx") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value after " + arg);
            }
            params.nx = std::stoul(argv[++i]);
        }
        else if (arg == "-ny" || arg == "--ny") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value after " + arg);
            }
            params.ny = std::stoul(argv[++i]);
        }
        else if (arg == "-nt" || arg == "--nt") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value after " + arg);
            }
            params.nt = std::stoul(argv[++i]);
        }
        else if (arg == "-dt" || arg == "--dt") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value after " + arg);
            }
            params.dt = std::stod(argv[++i]);
        }
        else if (arg == "-stabp" || arg == "--stabp") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value after " + arg);
            }
            params.stabp = std::stod(argv[++i]);
        }
        else if (arg == "-alpha" || arg == "--alpha") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value after " + arg);
            }
            params.alpha = std::stod(argv[++i]);
        }
        else if (arg == "-solver" || arg == "--solver") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value after " + arg);
            }
            params.solver_type = parse_solver_type(argv[++i]);
        }
        else if (arg == "-scheme" || arg == "--scheme") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value after " + arg);
            }
            params.scheme_type = parse_scheme_type(argv[++i]);
        }
        else if (arg == "-format" || arg == "--format") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value after " + arg);
            }
            params.output_format = parse_output_format(argv[++i]);
        }
        else if (arg == "-o" || arg == "--output") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing filename after " + arg);
            }
            params.output_prefix = argv[++i];
        }
        else if (arg == "-interval" || arg == "--interval") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value after " + arg);
            }
            params.output_interval = std::stoul(argv[++i]);
        }
        else if (arg == "-no-mpi" || arg == "--no-mpi") {
            params.use_mpi = false;
        }
        else if (!filename.empty() && arg[0] != '-') {
            // Treat as filename if not starting with '-'
            filename = arg;
        }
    }

    return filename;
}

void ParamReader::print_help(const std::string& program_name) {
    std::cout << "Usage: " << program_name << " [options] [param_file]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help              Show this help message\n";
    std::cout << "  -f, --file <file>       Parameter file (text or JSON format)\n";
    std::cout << "\nGrid parameters:\n";
    std::cout << "  -nx, --nx <value>       Number of grid points in x-direction\n";
    std::cout << "  -ny, --ny <value>       Number of grid points in y-direction\n";
    std::cout << "\nTime parameters:\n";
    std::cout << "  -nt, --nt <value>       Number of time steps\n";
    std::cout << "  -dt, --dt <value>       Time step size\n";
    std::cout << "  -stabp, --stabp <value> Stability parameter\n";
    std::cout << "\nPhysical parameters:\n";
    std::cout << "  -alpha, --alpha <value> Diffusion coefficient\n";
    std::cout << "\nSolver parameters:\n";
    std::cout << "  -solver, --solver <type> Solver type (Jacobi, SOR, CG, etc.)\n";
    std::cout << "  -scheme, --scheme <type> Time scheme (ImplicitEuler, CrankNicolson)\n";
    std::cout << "\nOutput parameters:\n";
    std::cout << "  -format, --format <fmt> Output format (Text, VTK, HDF5)\n";
    std::cout << "  -o, --output <prefix>   Output file prefix\n";
    std::cout << "  -interval, --interval N Output every N steps\n";
    std::cout << "\nMPI parameters:\n";
    std::cout << "  -no-mpi, --no-mpi       Disable MPI (run serial)\n";
    std::cout << "\nExample:\n";
    std::cout << "  " << program_name << " -f params.json\n";
    std::cout << "  " << program_name << " -nx 150 -ny 150 -nt 5 -stabp 0.1\n";
}

void ParamReader::print_params(const SimulationParams& params) {
    std::cout << "=== Simulation Parameters ===\n";
    std::cout << "Grid:\n";
    std::cout << "  nx = " << params.nx << "\n";
    std::cout << "  ny = " << params.ny << "\n";
    std::cout << "  lx = " << params.lx << "\n";
    std::cout << "  ly = " << params.ly << "\n";
    std::cout << "\nTime:\n";
    std::cout << "  dt = " << params.dt << "\n";
    std::cout << "  t_start = " << params.t_start << "\n";
    std::cout << "  t_end = " << params.t_end << "\n";
    std::cout << "  nt = " << params.nt << "\n";
    std::cout << "\nPhysical:\n";
    std::cout << "  alpha = " << params.alpha << "\n";
    std::cout << "  stabp = " << params.stabp << "\n";
    std::cout << "\nSolver:\n";
    std::cout << "  solver_type = " << params.get_solver_type_string() << "\n";
    std::cout << "  scheme_type = " << params.get_scheme_type_string() << "\n";
    std::cout << "  tolerance = " << params.solver_tolerance << "\n";
    std::cout << "  max_iterations = " << params.max_solver_iterations << "\n";
    std::cout << "\nOutput:\n";
    std::cout << "  format = " << params.get_output_format_string() << "\n";
    std::cout << "  prefix = " << params.output_prefix << "\n";
    std::cout << "  interval = " << params.output_interval << "\n";
    std::cout << "  output_initial = " << (params.output_initial ? "true" : "false") << "\n";
    std::cout << "  output_final = " << (params.output_final ? "true" : "false") << "\n";
    std::cout << "\nMPI:\n";
    std::cout << "  use_mpi = " << (params.use_mpi ? "true" : "false") << "\n";
    std::cout << "  num_procs_x = " << params.num_procs_x << "\n";
    std::cout << "  num_procs_y = " << params.num_procs_y << "\n";
    std::cout << "=============================\n";
}
