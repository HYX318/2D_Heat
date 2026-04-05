/**
 * @file solution_exporter.hpp
 * @brief Solution exporter for multiple output formats
 */

#ifndef SOLUTION_EXPORTER_HPP
#define SOLUTION_EXPORTER_HPP

#include <cstddef>
#include <string>
#include <fstream>
#include <stdexcept>
#include "../mesh/mesh2d.hpp"
#include "param_reader.hpp"

/**
 * @class SolutionExporter
 * @brief Multi-format solution exporter
 *
 * This class exports simulation solutions to various formats:
 * - Text format (backward compatible): I J x y u(x,y)
 * - VTK format: For visualization with ParaView, VisIt, etc.
 * - HDF5 format: For large datasets and scientific workflows
 *
 * Features:
 * - Automatic filename generation
 * - MPI-aware output (only root process writes global data)
 * - Time-series support
 * - Metadata preservation
 * - Error handling
 */
class SolutionExporter {
public:
    /**
     * @brief Constructor
     * @param output_prefix Output file prefix
     * @param output_format Output format
     * @param use_mpi Enable MPI mode (only root writes)
     */
    SolutionExporter(const std::string& output_prefix,
                      OutputFormat output_format,
                      bool use_mpi = false);

    /**
     * @brief Destructor
     */
    ~SolutionExporter();

    // Delete copy constructor and copy assignment
    SolutionExporter(const SolutionExporter&) = delete;
    SolutionExporter& operator=(const SolutionExporter&) = delete;

    // Move constructor and move assignment
    SolutionExporter(SolutionExporter&& other) noexcept = default;
    SolutionExporter& operator=(SolutionExporter&& other) noexcept = default;

    /**
     * @brief Export solution to file
     * @param mesh Mesh to export
     * @param time Current time
     * @param step Time step number
     * @param filename_suffix Optional filename suffix
     * @throws std::runtime_error if export fails
     */
    void export_solution(const Mesh2D& mesh,
                          double time,
                          size_t step,
                          const std::string& filename_suffix = "");

    /**
     * @brief Export solution to specific format
     * @param mesh Mesh to export
     * @param time Current time
     * @param step Time step number
     * @param format Output format to use
     * @param filename Output filename
     * @throws std::runtime_error if export fails
     */
    void export_solution(const Mesh2D& mesh,
                          double time,
                          size_t step,
                          OutputFormat format,
                          const std::string& filename);

    /**
     * @brief Export initial condition
     * @param mesh Mesh to export
     * @throws std::runtime_error if export fails
     */
    void export_initial(const Mesh2D& mesh);

    /**
     * @brief Export final solution
     * @param mesh Mesh to export
     * @param time Final time
     * @throws std::runtime_error if export fails
     */
    void export_final(const Mesh2D& mesh, double time);

    /**
     * @brief Set output prefix
     * @param prefix New output prefix
     */
    void set_output_prefix(const std::string& prefix) {
        output_prefix_ = prefix;
    }

    /**
     * @brief Get output prefix
     * @return Current output prefix
     */
    const std::string& get_output_prefix() const {
        return output_prefix_;
    }

    /**
     * @brief Set output format
     * @param format New output format
     */
    void set_output_format(OutputFormat format) {
        output_format_ = format;
    }

    /**
     * @brief Get output format
     * @return Current output format
     */
    OutputFormat get_output_format() const {
        return output_format_;
    }

    /**
     * @brief Generate filename for given step
     * @param step Time step number
     * @param suffix Optional suffix
     * @return Generated filename
     */
    std::string generate_filename(size_t step,
                                    const std::string& suffix = "") const;

    /**
     * @brief Get format-specific extension
     * @param format Output format
     * @return File extension (including dot)
     */
    static std::string get_format_extension(OutputFormat format);

    /**
     * @brief Get format name
     * @param format Output format
     * @return Format name as string
     */
    static std::string get_format_name(OutputFormat format);

private:
    std::string output_prefix_;  ///< Output file prefix
    OutputFormat output_format_; ///< Current output format
    bool use_mpi_;              ///< MPI mode flag

    /**
     * @brief Export to text format
     * @param mesh Mesh to export
     * @param time Current time
     * @param step Time step number
     * @param filename Output filename
     */
    void export_text(const Mesh2D& mesh,
                      double time,
                      size_t step,
                      const std::string& filename);

    /**
     * @brief Export to VTK format
     * @param mesh Mesh to export
     * @param time Current time
     * @param step Time step number
     * @param filename Output filename
     */
    void export_vtk(const Mesh2D& mesh,
                     double time,
                     size_t step,
                     const std::string& filename);

    /**
     * @brief Export to HDF5 format
     * @param mesh Mesh to export
     * @param time Current time
     * @param step Time step number
     * @param filename Output filename
     */
    void export_hdf5(const Mesh2D& mesh,
                       double time,
                       size_t step,
                       const std::string& filename);

    /**
     * @brief Check if file can be opened
     * @param filename File to check
     * @return true if file can be opened, false otherwise
     */
    bool can_open_file(const std::string& filename) const;

    /**
     * @brief Create directory for output if needed
     * @param filename Output filename
     */
    void ensure_directory_exists(const std::string& filename) const;
};

#endif // SOLUTION_EXPORTER_HPP
