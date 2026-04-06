/**
 * @file solution_exporter.cpp
 * @brief Implementation of solution exporter
 */

#include "solution_exporter.hpp"
#include <iomanip>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>

SolutionExporter::SolutionExporter(const std::string& output_prefix,
                                     OutputFormat output_format,
                                     bool use_mpi)
    : output_prefix_(output_prefix)
    , output_format_(output_format)
    , use_mpi_(use_mpi)
{
    // Logger not available - can uncomment when logging is properly set up
    // utils::Logger::get_instance().log_info(
    //     "SolutionExporter initialized with format: " + get_format_name(output_format_)
    // );
}

SolutionExporter::~SolutionExporter() = default;

void SolutionExporter::export_solution(const Mesh2D& mesh,
                                        double time,
                                        size_t step,
                                        const std::string& filename_suffix) {
    std::string filename = generate_filename(step, filename_suffix);
    export_solution(mesh, time, step, output_format_, filename);
}

void SolutionExporter::export_solution(const Mesh2D& mesh,
                                        double time,
                                        size_t step,
                                        OutputFormat format,
                                        const std::string& filename) {
    // Logger not available - can uncomment when logging is properly set up
    // utils::Logger::get_instance().log_info(
    //     "Exporting solution to " + filename + " (format: " + get_format_name(format) + ")"
    // );

    // Check if we should write (MPI mode: only root writes)
    if (use_mpi_) {
#ifdef USE_MPI
        // In MPI mode, check if we're root
        // This would require MPI context - for now, we'll write on all processes
        // In a real implementation, we'd gather data to root and write there
#else
        // No MPI, always write
#endif
    }

    // Ensure directory exists
    ensure_directory_exists(filename);

    // Export based on format
    switch (format) {
        case OutputFormat::Text:
            export_text(mesh, time, step, filename);
            break;
        case OutputFormat::VTK:
            export_vtk(mesh, time, step, filename);
            break;
        case OutputFormat::HDF5:
            export_hdf5(mesh, time, step, filename);
            break;
        default:
            throw std::runtime_error("Unknown output format");
    }
}

void SolutionExporter::export_initial(const Mesh2D& mesh) {
    export_solution(mesh, 0.0, 0, "initial");
}

void SolutionExporter::export_final(const Mesh2D& mesh, double time) {
    export_solution(mesh, time, std::numeric_limits<size_t>::max(), "final");
}

std::string SolutionExporter::generate_filename(size_t step,
                                                   const std::string& suffix) const {
    std::ostringstream oss;
    oss << output_prefix_;

    if (!suffix.empty()) {
        oss << "_" << suffix;
    }

    if (step != 0 && step != std::numeric_limits<size_t>::max()) {
        oss << "_step" << std::setw(6) << std::setfill('0') << step;
    }

    oss << get_format_extension(output_format_);
    return oss.str();
}

std::string SolutionExporter::get_format_extension(OutputFormat format) {
    switch (format) {
        case OutputFormat::Text: return ".txt";
        case OutputFormat::VTK: return ".vtk";
        case OutputFormat::HDF5: return ".h5";
        default: return ".out";
    }
}

std::string SolutionExporter::get_format_name(OutputFormat format) {
    switch (format) {
        case OutputFormat::Text: return "Text";
        case OutputFormat::VTK: return "VTK";
        case OutputFormat::HDF5: return "HDF5";
        default: return "Unknown";
    }
}

void SolutionExporter::export_text(const Mesh2D& mesh,
                                    double time,
                                    size_t step,
                                    const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    // Write header
    file << "# Solution at time = " << time << ", step = " << step << "\n";
    file << "# Grid: " << mesh.nx() << " x " << mesh.ny() << "\n";
    file << "# Physical domain: [0, " << mesh.lx() << "] x [0, " << mesh.ly() << "]\n";
    file << "# Format: I J x y u(x,y)\n";

    // Write data
    file << std::scientific << std::setprecision(10);
    for (size_t j = 0; j < mesh.ny(); ++j) {
        for (size_t i = 0; i < mesh.nx(); ++i) {
            double x = mesh.x_coord(i);
            double y = mesh.y_coord(j);
            double u = mesh(i, j);
            file << i << " " << j << " " << x << " " << y << " " << u << "\n";
        }
    }

    file.close();

    // Logger not available - can uncomment when logging is properly set up
    // utils::Logger::get_instance().log_info(
    //     "Exported " + std::to_string(mesh.nx() * mesh.ny()) +
    //     " points to " + filename
    // );
}

void SolutionExporter::export_vtk(const Mesh2D& mesh,
                                   double time,
                                   size_t step,
                                   const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    size_t nx = mesh.nx();
    size_t ny = mesh.ny();
    size_t num_points = nx * ny;

    // VTK header
    file << "# vtk DataFile Version 3.0\n";
    file << "Heat Equation Solution at time = " << time << ", step = " << step << "\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_GRID\n";

    // Grid dimensions
    file << "DIMENSIONS " << nx << " " << ny << " 1\n";

    // Points
    file << "POINTS " << num_points << " double\n";
    file << std::scientific << std::setprecision(10);
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            double x = mesh.x_coord(i);
            double y = mesh.y_coord(j);
            file << x << " " << y << " 0.0\n";
        }
    }

    // Point data
    file << "POINT_DATA " << num_points << "\n";

    // Scalar field (temperature)
    file << "SCALARS temperature double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            file << mesh(i, j) << "\n";
        }
    }

    file.close();

    // Logger not available - can uncomment when logging is properly set up
    // utils::Logger::get_instance().log_info(
    //     "Exported VTK file: " + filename +
    //     " " + std::to_string(num_points) + " points)"
    // );
}

void SolutionExporter::export_hdf5(const Mesh2D& mesh,
                                     double time,
                                     size_t step,
                                     const std::string& filename) {
    // Note: This is a placeholder for HDF5 export
    // In a real implementation, we would use the HDF5 C++ API
    // or high-level libraries like h5pp

    // Logger not available - can uncomment when logging is properly set up
    // utils::Logger::get_instance().log_warning(
    //     "HDF5 export not implemented - falling back to text format"
    // );

    // Fall back to text format
    std::string txt_filename = filename;
    if (txt_filename.size() > 3) {
        txt_filename = txt_filename.substr(0, txt_filename.size() - 3) + ".txt";
    }
    export_text(mesh, time, step, txt_filename);
}

bool SolutionExporter::can_open_file(const std::string& filename) const {
    std::ifstream file(filename);
    return file.good();
}

void SolutionExporter::ensure_directory_exists(const std::string& filename) const {
    // Find of last directory separator
    size_t pos = filename.find_last_of('/');
    if (pos == std::string::npos) {
        pos = filename.find_last_of('\\');
    }

    if (pos != std::string::npos) {
        std::string dir = filename.substr(0, pos);
        struct stat st;
        if (stat(dir.c_str(), &st) != 0) {
            // Directory doesn't exist, create it
            #ifdef _WIN32
            mkdir(dir.c_str());
            #else
            mkdir(dir.c_str(), 0755);
            #endif
        }
    }
}
