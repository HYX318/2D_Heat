/**
 * @file param_reader.hpp
 * @brief Parameter reader for simulation configuration
 */

#ifndef PARAM_READER_HPP
#define PARAM_READER_HPP

#include <cstddef>
#include <string>
#include <map>
#include <stdexcept>
#include <vector>
#include "../core/solver/solver_interface.hpp"
#include "../core/scheme/time_scheme_interface.hpp"

/**
 * @enum OutputFormat
 * @brief Supported output formats for solution files
 */
enum class OutputFormat {
    Text,   ///< Plain text format (backward compatible)
    VTK,    ///< VTK format for visualization
    HDF5    ///< HDF5 format for large datasets
};

/**
 * @enum InputFormat
 * @brief Supported input formats for parameter files
 */
enum class InputFormat {
    Auto,   ///< Auto-detect format
    Text,   ///< Plain text format
    JSON    ///< JSON format
};

/**
 * @struct SimulationParams
 * @brief All simulation parameters
 */
struct SimulationParams {
    // Grid parameters
    size_t nx;                      ///< Number of grid points in x-direction
    size_t ny;                      ///< Number of grid points in y-direction
    double lx;                      ///< Physical length in x-direction
    double ly;                      ///< Physical length in y-direction

    // Time parameters
    double dt;                      ///< Time step size
    double t_start;                 ///< Start time
    double t_end;                   ///< End time
    size_t nt;                      ///< Number of time steps

    // Physical parameters
    double alpha;                   ///< Diffusion coefficient
    double stabp;                   ///< Stability parameter (for validation)

    // Solver parameters
    SolverType solver_type;         ///< Type of linear solver
    TimeSchemeType scheme_type;     ///< Type of time integration scheme
    double solver_tolerance;        ///< Solver tolerance
    size_t max_solver_iterations;   ///< Maximum solver iterations

    // Output parameters
    OutputFormat output_format;     ///< Output file format
    std::string output_prefix;      ///< Output file prefix
    size_t output_interval;         ///< Output every N steps
    bool output_initial;            ///< Output initial condition
    bool output_final;              ///< Output final solution

    // MPI parameters
    bool use_mpi;                   ///< Enable MPI
    int num_procs_x;                ///< Number of processes in x-direction (0 = auto)
    int num_procs_y;                ///< Number of processes in y-direction (0 = auto)

    /**
     * @brief Default constructor with sensible defaults
     */
    SimulationParams()
        : nx(100)
        , ny(100)
        , lx(1.0)
        , ly(1.0)
        , dt(0.001)
        , t_start(0.0)
        , t_end(1.0)
        , nt(1000)
        , alpha(1.0)
        , stabp(0.25)
        , solver_type(SolverType::Jacobi)
        , scheme_type(TimeSchemeType::ImplicitEuler)
        , solver_tolerance(1e-6)
        , max_solver_iterations(10000)
        , output_format(OutputFormat::Text)
        , output_prefix("solution")
        , output_interval(100)
        , output_initial(true)
        , output_final(true)
        , use_mpi(true)
        , num_procs_x(0)
        , num_procs_y(0)
    {}

    /**
     * @brief Validate parameters
     * @throws std::invalid_argument if parameters are invalid
     */
    void validate() const;

    /**
     * @brief Compute time step from stability parameter
     * @return Time step size
     */
    double compute_dt() const;

    /**
     * @brief Compute stability parameter from time step
     * @return Stability parameter
     */
    double compute_stabp() const;

    /**
     * @brief Get solver type string
     * @return String representation
     */
    std::string get_solver_type_string() const;

    /**
     * @brief Get scheme type string
     * @return String representation
     */
    std::string get_scheme_type_string() const;

    /**
     * @brief Get output format string
     * @return String representation
     */
    std::string get_output_format_string() const;
};

/**
 * @class ParamReader
 * @brief Parameter reader for simulation configuration
 *
 * This class reads simulation parameters from various input formats:
 * - Text format (backward compatible with old Param.in format)
 * - JSON format (new, more structured)
 *
 * Features:
 * - Automatic format detection
 * - Parameter validation
 * - Default value support
 * - Error handling with descriptive messages
 * - Command line argument support
 */
class ParamReader {
public:
    /**
     * @brief Constructor
     */
    ParamReader();

    /**
     * @brief Destructor
     */
    ~ParamReader() = default;

    /**
     * @brief Read parameters from file
     * @param filename Input file name
     * @param format Input format (Auto to auto-detect)
     * @return Simulation parameters
     * @throws std::runtime_error if file cannot be read
     * @throws std::invalid_argument if parameters are invalid
     */
    SimulationParams read_from_file(const std::string& filename,
                                     InputFormat format = InputFormat::Auto);

    /**
     * @brief Read parameters from text format file
     * @param filename Input file name
     * @return Simulation parameters
     * @throws std::runtime_error if file cannot be read
     */
    SimulationParams read_text_format(const std::string& filename);

    /**
     * @brief Read parameters from JSON format file
     * @param filename Input file name
     * @return Simulation parameters
     * @throws std::runtime_error if file cannot be read or JSON is invalid
     */
    SimulationParams read_json_format(const std::string& filename);

    /**
     * @brief Parse command line arguments
     * @param argc Argument count
     * @param argv Argument vector
     * @param params Parameters to update (updated in place)
     * @return Filename if provided, empty string otherwise
     * @throws std::invalid_argument if arguments are invalid
     */
    std::string parse_command_line(int argc, char** argv,
                                    SimulationParams& params);

    /**
     * @brief Print help message
     * @param program_name Program name
     */
    static void print_help(const std::string& program_name);

    /**
     * @brief Print parameters
     * @param params Parameters to print
     */
    static void print_params(const SimulationParams& params);

private:
    /**
     * @brief Detect file format from extension
     * @param filename File name
     * @return Detected format
     */
    InputFormat detect_format(const std::string& filename) const;

    /**
     * @brief Parse solver type from string
     * @param str String representation
     * @return Solver type
     * @throws std::invalid_argument if string is invalid
     */
    SolverType parse_solver_type(const std::string& str) const;

    /**
     * @brief Parse scheme type from string
     * @param str String representation
     * @return Scheme type
     * @throws std::invalid_argumentation if string is invalid
     */
    TimeSchemeType parse_scheme_type(const std::string& str) const;

    /**
     * @brief Parse output format from string
     * @param str String representation
     * @return Output format
     * @throws std::invalid_argument if string is invalid
     */
    OutputFormat parse_output_format(const std::string& str) const;

    /**
     * @brief Trim whitespace from string
     * @param str Input string
     * @return Trimmed string
     */
    std::string trim(const std::string& str) const;

    /**
     * @brief Split string by delimiter
     * @param str Input string
     * @param delim Delimiter
     * @return Vector of tokens
     */
    std::vector<std::string> split(const std::string& str,
                                     char delim) const;
};

#endif // PARAM_READER_HPP
