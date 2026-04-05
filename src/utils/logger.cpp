#include "logger.hpp"
#include <ctime>

// Constructor for console output
Logger::Logger(LogLevel min_level)
    : min_level_(min_level)
    , file_stream_(nullptr)
    , timestamps_enabled_(true)
    , mpi_rank_enabled_(false)
    , mpi_rank_(-1)
    , use_color_(true)
{
}

// Constructor for file output
Logger::Logger(const std::string& filename, LogLevel min_level)
    : min_level_(min_level)
    , file_stream_(new std::ofstream(filename))
    , timestamps_enabled_(true)
    , mpi_rank_enabled_(false)
    , mpi_rank_(-1)
    , use_color_(false)  // No color for file output
{
    if (!file_stream_->is_open()) {
        std::cerr << "Failed to open log file: " << filename << std::endl;
        file_stream_.reset();
    }
}

// Destructor
Logger::~Logger()
{
    if (file_stream_ && file_stream_->is_open()) {
        file_stream_->close();
    }
}

// Check if a log level should be displayed
bool Logger::should_log(LogLevel level) const
{
    return level >= min_level_;
}

// Convert log level to string
std::string Logger::level_to_string(LogLevel level) const
{
    switch (level) {
        case LogLevel::DEBUG:
            return "DEBUG";
        case LogLevel::INFO:
            return "INFO";
        case LogLevel::WARNING:
            return "WARNING";
        case LogLevel::ERROR:
            return "ERROR";
        case LogLevel::FATAL:
            return "FATAL";
        default:
            return "UNKNOWN";
    }
}

// Get ANSI color code for log level
std::string Logger::level_to_color(LogLevel level) const
{
    if (!use_color_) {
        return "";
    }

    switch (level) {
        case LogLevel::DEBUG:
            return "\x1b[36m";    // Cyan
        case LogLevel::INFO:
            return "\x1b[32m";     // Green
        case LogLevel::WARNING:
            return "\x1b[33m";     // Yellow
        case LogLevel::ERROR:
            return "\x1b[31m";     // Red
        case LogLevel::FATAL:
            return "\x1b[35m";     // Magenta
        default:
            return "\x1b[0m";      // Reset
    }
}

// Get current timestamp as formatted string
std::string Logger::get_timestamp() const
{
    auto now = std::chrono::system_clock::now();
    auto time_t_now = std::chrono::system_clock::to_time_t(now);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        now.time_since_epoch()
    ) % 1000;

    std::tm* tm_time = std::localtime(&time_t_now);

    std::ostringstream ss;
    ss << std::setfill('0') << std::setw(4) << (tm_time->tm_year + 1900) << "-";
    ss << std::setfill('0') << std::setw(2) << (tm_time->tm_mon + 1) << "-";
    ss << std::setfill('0') << std::setw(2) << tm_time->tm_mday << " ";
    ss << std::setfill('0') << std::setw(2) << tm_time->tm_hour << ":";
    ss << std::setfill('0') << std::setw(2) << tm_time->tm_min << ":";
    ss << std::setfill('0') << std::setw(2) << tm_time->tm_sec << ".";
    ss << std::setfill('0') << std::setw(3) << ms.count();

    return ss.str();
}

// Format a log message
std::string Logger::format_message(LogLevel level, const std::string& category, const std::string& message) const
{
    std::ostringstream ss;

    // Add timestamp if enabled
    if (timestamps_enabled_) {
        ss << "[" << get_timestamp() << "] ";
    }

    // Add MPI rank if enabled
    if (mpi_rank_enabled_ && mpi_rank_ >= 0) {
        ss << "[RANK-" << mpi_rank_ << "] ";
    }

    // Add level
    ss << "[" << level_to_string(level) << "] ";

    // Add category if provided
    if (!category.empty()) {
        ss << "[" << category << "] ";
    }

    // Add message
    ss << message;

    return ss.str();
}

// Write message to output stream
void Logger::write_message(const std::string& formatted_message, LogLevel level)
{
    std::lock_guard<std::mutex> lock(mutex_);

    std::string reset_color = use_color_ ? "\x1b[0m" : "";
    std::string color_code = level_to_color(level);

    if (file_stream_ && file_stream_->is_open()) {
        // Write to file (no colors)
        *file_stream_ << formatted_message << std::endl;
        file_stream_->flush();
    } else {
        // Write to console (with colors)
        std::cout << color_code << formatted_message << reset_color << std::endl;
    }
}

// Log a message with specified level
void Logger::log(LogLevel level, const std::string& message)
{
    log(level, "", message);
}

// Log a message with specified level and category
void Logger::log(LogLevel level, const std::string& category, const std::string& message)
{
    if (!should_log(level)) {
        return;
    }

    std::string formatted_message = format_message(level, category, message);
    write_message(formatted_message, level);
}

// DEBUG level log
void Logger::debug(const std::string& message)
{
    log(LogLevel::DEBUG, message);
}

// INFO level log
void Logger::info(const std::string& message)
{
    log(LogLevel::INFO, message);
}

// WARNING level log
void Logger::warning(const std::string& message)
{
    log(LogLevel::WARNING, message);
}

// ERROR level log
void Logger::error(const std::string& message)
{
    log(LogLevel::ERROR, message);
}

// FATAL level log
void Logger::fatal(const std::string& message)
{
    log(LogLevel::FATAL, message);
}

// Set minimum log level
void Logger::set_level(LogLevel level)
{
    std::lock_guard<std::mutex> lock(mutex_);
    min_level_ = level;
}

// Enable or disable timestamps
void Logger::enable_timestamps(bool enable)
{
    std::lock_guard<std::mutex> lock(mutex_);
    timestamps_enabled_ = enable;
}

// Enable or disable MPI rank output
void Logger::enable_mpi_rank(bool enable, int rank)
{
    std::lock_guard<std::mutex> lock(mutex_);
    mpi_rank_enabled_ = enable;
    mpi_rank_ = rank;
}
