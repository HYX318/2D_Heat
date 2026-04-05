#include "logger.hpp"
#include <thread>
#include <vector>
#include <chrono>

void test_logger_console() {
    std::cout << "=== Testing Console Logger ===" << std::endl;

    Logger logger(Logger::LogLevel::DEBUG);
    logger.enable_mpi_rank(true, 0);
    logger.set_level(Logger::LogLevel::INFO);

    logger.info("Simulation started");
    logger.debug("This debug message should not appear (level too low)");
    logger.warning("Convergence slower than expected");
    logger.error("Failed to read parameter file");
    logger.fatal("Critical error, aborting");

    logger.log(Logger::LogLevel::INFO, "IO", "File loaded successfully");
    logger.log(Logger::LogLevel::WARNING, "Solver", "Using default parameters");

    std::cout << std::endl;
}

void test_logger_file() {
    std::cout << "=== Testing File Logger ===" << std::endl;

    Logger logger("test_log.txt", Logger::LogLevel::DEBUG);
    logger.enable_mpi_rank(true, 1);
    logger.enable_timestamps(true);

    logger.debug("Debug message to file");
    logger.info("Info message to file");
    logger.warning("Warning message to file");
    logger.error("Error message to file");
    logger.fatal("Fatal message to file");

    logger.log(Logger::LogLevel::INFO, "MPI", "Rank 1 initialized");

    std::cout << "Messages written to test_log.txt" << std::endl << std::endl;
}

void test_logger_thread_safety() {
    std::cout << "=== Testing Thread Safety ===" << std::endl;

    Logger logger(Logger::LogLevel::DEBUG);
    logger.set_level(Logger::LogLevel::INFO);

    std::vector<std::thread> threads;
    for (int i = 0; i < 5; ++i) {
        threads.emplace_back([&logger, i]() {
            for (int j = 0; j < 10; ++j) {
                logger.info("Thread " + std::to_string(i) + " message " + std::to_string(j));
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }
        });
    }

    for (auto& thread : threads) {
        thread.join();
    }

    std::cout << "Thread safety test completed" << std::endl << std::endl;
}

void test_logger_level_filtering() {
    std::cout << "=== Testing Level Filtering ===" << std::endl;

    Logger logger(Logger::LogLevel::INFO);

    std::cout << "With level set to INFO:" << std::endl;
    logger.set_level(Logger::LogLevel::INFO);
    logger.debug("Debug (hidden)");
    logger.info("Info (visible)");
    logger.warning("Warning (visible)");
    logger.error("Error (visible)");
    logger.fatal("Fatal (visible)");

    std::cout << "\nWith level set to ERROR:" << std::endl;
    logger.set_level(Logger::LogLevel::ERROR);
    logger.debug("Debug (hidden)");
    logger.info("Info (hidden)");
    logger.warning("Warning (hidden)");
    logger.error("Error (visible)");
    logger.fatal("Fatal (visible)");

    std::cout << "\nWith level set to DEBUG:" << std::endl;
    logger.set_level(Logger::LogLevel::DEBUG);
    logger.debug("Debug (visible)");
    logger.info("Info (visible)");
    logger.warning("Warning (visible)");
    logger.error("Error (visible)");
    logger.fatal("Fatal (visible)");

    std::cout << std::endl;
}

void test_logger_features() {
    std::cout << "=== Testing Features ===" << std::endl;

    Logger logger(Logger::LogLevel::INFO);

    std::cout << "With timestamps:" << std::endl;
    logger.enable_timestamps(true);
    logger.enable_mpi_rank(false);
    logger.info("Message with timestamp");

    std::cout << "\nWith MPI rank:" << std::endl;
    logger.enable_timestamps(false);
    logger.enable_mpi_rank(true, 42);
    logger.info("Message with rank");

    std::cout << "\nWith both:" << std::endl;
    logger.enable_timestamps(true);
    logger.enable_mpi_rank(true, 7);
    logger.info("Message with both");

    std::cout << "\nWith category:" << std::endl;
    logger.log(Logger::LogLevel::INFO, "Network", "Connection established");

    std::cout << std::endl;
}

int main() {
    test_logger_console();
    test_logger_file();
    test_logger_thread_safety();
    test_logger_level_filtering();
    test_logger_features();

    std::cout << "All tests completed successfully!" << std::endl;
    return 0;
}
