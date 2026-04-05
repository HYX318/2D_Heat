#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <mutex>
#include <memory>
#include <chrono>
#include <iomanip>
#include <sstream>

/**
 * @brief Logger class for structured logging with different levels and output targets
 *
 * Features:
 * - Thread-safe logging with mutex protection
 * - Multiple log levels (DEBUG, INFO, WARNING, ERROR, FATAL)
 * - Timestamp support
 * - MPI rank support for parallel debugging
 * - Multiple output targets (console, file)
 * - Colored console output
 */
class Logger {
public:
    /**
     * @brief Log level enumeration
     */
    enum class LogLevel {
        DEBUG,
        INFO,
        WARNING,
        ERROR,
        FATAL
    };

    /**
     * @brief Construct a Logger that outputs to console
     * @param min_level Minimum log level to be displayed
     */
    explicit Logger(LogLevel min_level = LogLevel::INFO);

    /**
     * @brief Construct a Logger that outputs to a file
     * @param filename Output file path
     * @param min_level Minimum log level to be displayed
     */
    Logger(const std::string& filename, LogLevel min_level = LogLevel::INFO);

    /**
     * @brief Destructor
     */
    ~Logger();

    // Delete copy constructor and copy assignment for thread safety
    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;

    // Delete move constructor and move assignment (mutex is not movable)
    Logger(Logger&&) = delete;
    Logger& operator=(Logger&&) = delete;

    /**
     * @brief Log a DEBUG level message
     * @param message Message to log
     */
    void debug(const std::string& message);

    /**
     * @brief Log an INFO level message
     * @param message Message to log
     */
    void info(const std::string& message);

    /**
     * @brief Log a WARNING level message
     * @param message Message to log
     */
    void warning(const std::string& message);

    /**
     * @brief Log an ERROR level message
     * @param message Message to log
     */
    void error(const std::string& message);

    /**
     * @brief Log a FATAL level message
     * @param message Message to log
     */
    void fatal(const std::string& message);

    /**
     * @brief Log a message with specified level and category
     * @param level Log level
     * @param category Category identifier
     * @param message Message to log
     */
    void log(LogLevel level, const std::string& category, const std::string& message);

    /**
     * @brief Log a message with specified level
     * @param level Log level
     * @param message Message to log
     */
    void log(LogLevel level, const std::string& message);

    /**
     * @brief Set the minimum log level
     * @param level Minimum level to be displayed
     */
    void set_level(LogLevel level);

    /**
     * @brief Enable or disable timestamps
     * @param enable true to enable timestamps
     */
    void enable_timestamps(bool enable);

    /**
     * @brief Enable or disable MPI rank output
     * @param enable true to enable MPI rank output
     * @param rank MPI rank (-1 for no rank)
     */
    void enable_mpi_rank(bool enable, int rank = -1);

private:
    /**
     * @brief Get the string representation of a log level
     * @param level Log level
     * @return String representation
     */
    std::string level_to_string(LogLevel level) const;

    /**
     * @brief Get the color code for a log level (for console output)
     * @param level Log level
     * @return ANSI color code
     */
    std::string level_to_color(LogLevel level) const;

    /**
     * @brief Get the current timestamp as a formatted string
     * @return Formatted timestamp string
     */
    std::string get_timestamp() const;

    /**
     * @brief Format a log message
     * @param level Log level
     * @param category Category (optional)
     * @param message Message content
     * @return Formatted message string
     */
    std::string format_message(LogLevel level, const std::string& category, const std::string& message) const;

    /**
     * @brief Write a message to the output stream
     * @param formatted_message Formatted message
     * @param level Log level (for color coding)
     */
    void write_message(const std::string& formatted_message, LogLevel level);

    /**
     * @brief Check if a log level should be displayed
     * @param level Log level to check
     * @return true if level should be displayed
     */
    bool should_log(LogLevel level) const;

    // Member variables
    LogLevel min_level_;
    std::unique_ptr<std::ofstream> file_stream_;
    std::mutex mutex_;
    bool timestamps_enabled_;
    bool mpi_rank_enabled_;
    int mpi_rank_;
    bool use_color_;
};

#endif // LOGGER_HPP
