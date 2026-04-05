#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>
#include <string>
#include <iostream>

/**
 * @brief High-resolution timer for performance measurement
 *
 * This class provides precise timing functionality with support for:
 * - Named timers
 * - Lap timing
 * - Hierarchical timing
 * - RAII-based automatic timing (via ScopedTimer)
 *
 * Uses std::chrono::high_resolution_clock for maximum precision.
 */
class Timer {
public:
    using clock_type = std::chrono::high_resolution_clock;
    using time_point = clock_type::time_point;
    using duration = clock_type::duration;

    /**
     * @brief Construct an unnamed timer and immediately start timing
     */
    Timer()
        : name_("Timer")
        , start_time_(clock_type::now())
        , stop_time_(time_point{})
        , lap_count_(0)
        , last_lap_time_(start_time_)
        , running_(true)
    {}

    /**
     * @brief Construct a named timer and immediately start timing
     * @param name Name of the timer
     */
    explicit Timer(const std::string& name)
        : name_(name)
        , start_time_(clock_type::now())
        , stop_time_(time_point{})
        , lap_count_(0)
        , last_lap_time_(start_time_)
        , running_(true)
    {}

    /**
     * @brief Start or restart timing
     *
     * If the timer is already running, this resets the start time.
     */
    void start() {
        start_time_ = clock_type::now();
        last_lap_time_ = start_time_;
        stop_time_ = time_point{};
        running_ = true;
    }

    /**
     * @brief Stop timing
     *
     * Stores the stop time for later elapsed time calculations.
     */
    void stop() {
        stop_time_ = clock_type::now();
        running_ = false;
    }

    /**
     * @brief Reset the timer to initial state
     *
     * Resets all timing information and lap counts.
     * Does not start the timer automatically.
     */
    void reset() {
        start_time_ = time_point{};
        stop_time_ = time_point{};
        lap_count_ = 0;
        last_lap_time_ = time_point{};
        running_ = false;
    }

    /**
     * @brief Check if the timer is currently running
     * @return true if running, false otherwise
     */
    bool is_running() const {
        return running_;
    }

    /**
     * @brief Get elapsed time in seconds
     * @return Elapsed time in seconds as a double
     */
    double elapsed() const {
        if (running_) {
            return std::chrono::duration<double>(clock_type::now() - start_time_).count();
        } else {
            return std::chrono::duration<double>(stop_time_ - start_time_).count();
        }
    }

    /**
     * @brief Get elapsed time in milliseconds
     * @return Elapsed time in milliseconds as a double
     */
    double elapsed_ms() const {
        return elapsed() * 1000.0;
    }

    /**
     * @brief Get elapsed time in microseconds
     * @return Elapsed time in microseconds as a double
     */
    double elapsed_us() const {
        return elapsed() * 1000000.0;
    }

    /**
     * @brief Set the name of the timer
     * @param name New name for the timer
     */
    void set_name(const std::string& name) {
        name_ = name;
    }

    /**
     * @brief Get the name of the timer
     * @return The timer's name
     */
    const std::string& name() const {
        return name_;
    }

    /**
     * @brief Get the number of laps recorded
     * @return The lap count
     */
    size_t lap_count() const {
        return lap_count_;
    }

    /**
     * @brief Record a lap point and return the time since the last lap
     *
     * If the timer is not running, it will be started automatically.
     * The first lap returns the time from the start to the first lap call.
     * Subsequent laps return the time since the previous lap call.
     *
     * @return Time in seconds since the last lap (or start, for first lap)
     */
    double lap() {
        if (!running_) {
            start();
        }

        time_point current_time = clock_type::now();
        double lap_time = std::chrono::duration<double>(current_time - last_lap_time_).count();
        last_lap_time_ = current_time;
        ++lap_count_;

        return lap_time;
    }

    /**
     * @brief Reset the lap count
     *
     * Does not affect the timing information, only resets the lap counter.
     */
    void reset_lap_count() {
        lap_count_ = 0;
        last_lap_time_ = running_ ? clock_type::now() : time_point{};
    }

private:
    std::string name_;          ///< Timer name
    time_point start_time_;     ///< Time when timer was started
    time_point stop_time_;      ///< Time when timer was stopped
    size_t lap_count_;          ///< Number of laps recorded
    time_point last_lap_time_;  ///< Time of the last lap
    bool running_;              ///< Whether the timer is currently running
};

/**
 * @brief RAII-based scoped timer
 *
 * Creates a timer with the given name and automatically stops it
 * when it goes out of scope. Useful for automatic timing of code blocks.
 *
 * Example:
 * @code
 * {
 *     ScopedTimer timer("Computation");
 *     // ... do work ...
 * } // timer automatically stops here and prints elapsed time
 * @endcode
 */
class ScopedTimer {
public:
    /**
     * @brief Construct a scoped timer with a name
     * @param name Name of the timer
     * @param auto_print Whether to automatically print elapsed time on destruction
     */
    explicit ScopedTimer(const std::string& name, bool auto_print = true)
        : timer_(name)
        , auto_print_(auto_print)
    {}

    /**
     * @brief Construct a scoped timer that uses an existing timer
     * @param timer Reference to an existing timer
     * @param auto_print Whether to automatically print elapsed time on destruction
     */
    explicit ScopedTimer(Timer& timer, bool auto_print = true)
        : timer_(timer)
        , owns_timer_(false)
        , auto_print_(auto_print)
    {
        timer_.start();
    }

    /**
     * @brief Destroy the scoped timer
     *
     * Stops the timer and optionally prints the elapsed time.
     */
    ~ScopedTimer() {
        timer_.stop();
        if (auto_print_) {
            std::cout << "[" << timer_.name() << "] "
                      << timer_.elapsed() << " s ("
                      << timer_.elapsed_ms() << " ms)" << std::endl;
        }
    }

    // Delete copy constructor and copy assignment
    ScopedTimer(const ScopedTimer&) = delete;
    ScopedTimer& operator=(const ScopedTimer&) = delete;

    // Allow move operations
    ScopedTimer(ScopedTimer&& other) noexcept
        : timer_(std::move(other.timer_))
        , owns_timer_(other.owns_timer_)
        , auto_print_(other.auto_print_)
    {
        other.owns_timer_ = false;
        other.auto_print_ = false;
    }

    ScopedTimer& operator=(ScopedTimer&& other) noexcept {
        if (this != &other) {
            timer_ = std::move(other.timer_);
            owns_timer_ = other.owns_timer_;
            auto_print_ = other.auto_print_;
            other.owns_timer_ = false;
            other.auto_print_ = false;
        }
        return *this;
    }

    /**
     * @brief Get the underlying timer
     * @return Reference to the timer
     */
    Timer& timer() {
        return timer_;
    }

    /**
     * @brief Get the elapsed time in seconds
     * @return Elapsed time in seconds
     */
    double elapsed() const {
        return timer_.elapsed();
    }

    /**
     * @brief Get the elapsed time in milliseconds
     * @return Elapsed time in milliseconds
     */
    double elapsed_ms() const {
        return timer_.elapsed_ms();
    }

    /**
     * @brief Get the elapsed time in microseconds
     * @return Elapsed time in microseconds
     */
    double elapsed_us() const {
        return timer_.elapsed_us();
    }

private:
    Timer timer_;       ///< The timer object
    bool owns_timer_ = true;  ///< Whether this scoped timer owns the timer
    bool auto_print_;  ///< Whether to print on destruction
};

#endif // TIMER_HPP
