/**
 * @file timer_example.cpp
 * @brief Example demonstrating the usage of the Timer class
 */

#include "timer.hpp"
#include <iostream>
#include <thread>
#include <vector>

void simulate_work(int milliseconds) {
    std::this_thread::sleep_for(std::chrono::milliseconds(milliseconds));
}

int main() {
    std::cout << "=== Timer Class Examples ===" << std::endl << std::endl;

    // Example 1: Basic timing
    std::cout << "Example 1: Basic timing" << std::endl;
    Timer total_timer("Total Time");

    Timer step_timer("Step Timer");
    for (int i = 0; i < 3; ++i) {
        step_timer.start();
        simulate_work(100);
        step_timer.stop();
        std::cout << "  Step " << i << ": " << step_timer.elapsed_ms() << " ms" << std::endl;
    }

    total_timer.stop();
    std::cout << "  Total time: " << total_timer.elapsed() << " seconds" << std::endl;
    std::cout << "  Total time: " << total_timer.elapsed_ms() << " ms" << std::endl;
    std::cout << "  Total time: " << total_timer.elapsed_us() << " us" << std::endl;
    std::cout << std::endl;

    // Example 2: Lap timing
    std::cout << "Example 2: Lap timing" << std::endl;
    Timer loop_timer("Loop Timer");
    for (int i = 0; i < 3; ++i) {
        simulate_work(50);
        double lap_time = loop_timer.lap();
        std::cout << "  Lap " << i << ": " << lap_time << " s (" << lap_time * 1000 << " ms)" << std::endl;
    }
    std::cout << "  Total laps: " << loop_timer.lap_count() << std::endl;
    std::cout << std::endl;

    // Example 3: ScopedTimer with automatic printing
    std::cout << "Example 3: ScopedTimer (automatic)" << std::endl;
    {
        ScopedTimer timer("Scoped Block 1");
        simulate_work(200);
    } // timer automatically stops here and prints elapsed time
    std::cout << std::endl;

    // Example 4: ScopedTimer without auto-print
    std::cout << "Example 4: ScopedTimer (manual print)" << std::endl;
    {
        ScopedTimer timer("Scoped Block 2", false);
        simulate_work(150);
        std::cout << "  Inside scope: " << timer.elapsed_ms() << " ms" << std::endl;
    }
    std::cout << std::endl;

    // Example 5: Nested timers
    std::cout << "Example 5: Nested timing" << std::endl;
    Timer outer_timer("Outer");
    outer_timer.start();

    simulate_work(50);

    {
        Timer inner_timer("Inner");
        inner_timer.start();
        simulate_work(100);
        inner_timer.stop();
        std::cout << "  Inner time: " << inner_timer.elapsed_ms() << " ms" << std::endl;
    }

    outer_timer.stop();
    std::cout << "  Outer time: " << outer_timer.elapsed_ms() << " ms" << std::endl;
    std::cout << std::endl;

    // Example 6: Timer state operations
    std::cout << "Example 6: Timer state operations" << std::endl;
    Timer state_timer("State Test");
    std::cout << "  Initially running: " << (state_timer.is_running() ? "yes" : "no") << std::endl;
    std::cout << "  Initial name: " << state_timer.name() << std::endl;

    state_timer.stop();
    std::cout << "  After stop, running: " << (state_timer.is_running() ? "yes" : "no") << std::endl;

    state_timer.reset();
    std::cout << "  After reset, running: " << (state_timer.is_running() ? "yes" : "no") << std::endl;
    std::cout << "  After reset, elapsed: " << state_timer.elapsed() << " s" << std::endl;

    state_timer.set_name("Renamed Timer");
    std::cout << "  Renamed to: " << state_timer.name() << std::endl;
    std::cout << std::endl;

    // Example 7: Performance comparison
    std::cout << "Example 7: Performance measurement" << std::endl;
    const int n = 10000;

    Timer vector_timer("Vector allocation");
    vector_timer.start();
    std::vector<int> vec(n);
    for (int i = 0; i < n; ++i) {
        vec[i] = i * i;
    }
    vector_timer.stop();
    std::cout << "  Vector operations: " << vector_timer.elapsed_us() << " us" << std::endl;

    Timer array_timer("Array allocation");
    array_timer.start();
    int* arr = new int[n];
    for (int i = 0; i < n; ++i) {
        arr[i] = i * i;
    }
    delete[] arr;
    array_timer.stop();
    std::cout << "  Array operations: " << array_timer.elapsed_us() << " us" << std::endl;
    std::cout << std::endl;

    // Example 8: Continuous lap timing
    std::cout << "Example 8: Continuous lap timing with reset" << std::endl;
    Timer continuous_timer("Continuous Timer");
    for (int i = 0; i < 5; ++i) {
        simulate_work(30);
        double lap_time = continuous_timer.lap();
        std::cout << "  Lap " << i << ": " << lap_time * 1000 << " ms" << std::endl;
    }
    std::cout << "  Total laps: " << continuous_timer.lap_count() << std::endl;
    std::cout << "  Total time: " << continuous_timer.elapsed_ms() << " ms" << std::endl;

    continuous_timer.reset_lap_count();
    std::cout << "  After reset, lap count: " << continuous_timer.lap_count() << std::endl;

    // Continue timing after lap reset
    for (int i = 0; i < 3; ++i) {
        simulate_work(40);
        continuous_timer.lap();
    }
    std::cout << "  After more laps: " << continuous_timer.lap_count() << std::endl;
    std::cout << std::endl;

    std::cout << "=== All examples completed ===" << std::endl;

    return 0;
}
