# MPI Performance Profiling Guide

## Overview

This document describes the MPI performance profiling tools and integration available in the 2D Heat Equation project.

## Components

### 1. Profiler Class (`src/mpi/profiler.hpp`, `src/mpi/profiler.cpp`)

A comprehensive profiling tool for MPI applications that provides:

- **Hierarchical Timing**: Track total, communication, computation, and idle times
- **Communication Tracking**: Monitor message counts, data volumes, and operation types
- **Performance Metrics**: Calculate communication/computation ratios and load balance efficiency
- **Multiple Report Formats**: Console output, CSV, and JSON

#### Key Features

- Zero-overhead when disabled (compile-time option)
- Thread-safe statistics
- MPI Profiling Interface (PMPI) support
- Integration with MPI_Wtime() for accurate timing

### 2. Benchmark Suite (`tests/performance/`)

Comprehensive benchmarks for measuring MPI communication performance:

- **benchmark_mpi_communication.cpp**: Full MPI communication benchmark suite
  - Ghost cell exchange (blocking, async, sendrecv)
  - Reduction operations (reduce, allreduce)
  - Broadcast performance
  - All-to-all communication

### 3. Profiling Script (`scripts/profile_mpi.sh`)

Automated profiling script that supports:

- Built-in timer profiling
- Score-P instrumentation
- Intel VTune Amplifier
- Automatic report generation

## Usage

### Basic Profiler Usage

```cpp
#include "mpi/profiler.hpp"
#include "mpi/mpi_context.hpp"

int main(int argc, char** argv) {
    MPIContext mpi(argc, argv);

    // Create profiler (rank and size for reporting)
    Profiler profiler(true, mpi.rank(), mpi.size());

    // Profile computation
    profiler.start_computation("compute_step");
    // ... do computation ...
    profiler.end_computation("compute_step");

    // Profile communication
    profiler.start_communication("exchange_ghost");
    MPI_Isend(...);
    MPI_Irecv(...);
    MPI_Waitall(...);
    profiler.end_communication("exchange_ghost");

    // Record communication statistics
    profiler.record_send(dest, bytes, "ghost_exchange");
    profiler = profiler.record_recv(src, bytes, "ghost_exchange");

    // Print report
    profiler.print_report();

    // // Write reports
    profiler.write_report("report.txt");
    profiler.write_csv("data.csv");
    profiler.write_json("data.json");

    return 0;
}
```

### Running Benchmarks

### Basic benchmark execution:

```bash
# Build the project
mkdir build && cd build
cmake ..
make

# Run MPI communication benchmark
mpirun -np 4 bin/benchmark_mpi_communication

# Specify message sizes and iterations
mpirun -np 4 bin/benchmark_mpi_communication --sizes "1,10,100,1000,10000" --iters 100

# Specify output prefix
mpirun -np 4 bin/benchmark_mpi_communication --output my_results
```

### Using the Profiling Script

```bash
# Basic usage with built-in profiling
./scripts/profile_mpi.sh -p bin/benchmark_mpi_communication -n 4 -o results/ -v

# Use Score-P for detailed profiling
./scripts/profile_mpi.sh -p bin/heat_equation -n 4 -t scorep -o scorep_results/

# Use Intel VTune (if available)
./scripts/profile_mpi.sh -p bin/heat_equation -n 4 -t vtune -o vtune_results/

# With custom iterations and warmup
./scripts/profile_mpi.sh -p bin/benchmark_mpi_communication -n 8 -w 5 -m 10 -o results/
```

### Profiling Script Options

```
-n NUM_PROCES     Number of MPI processes (default: 4)
-o OUTPUT_DIR     Output directory for results (default: results/)
-p PROGRAM        MPI program to profile (required)
-a "ARGS"         Arguments to pass to the program
-t TOOL           Profiling tool: none, scorep, vtune (default: none)
-w WARMUP         Number of warmup runs (default: 1)
-m MEASURED       Number of measured runs (default: 3)
-v                Verbose output
-h                Show help message
```

## Report Formats

### Console Output

Human-readable summary printed to stdout:

```
===========================================
MPI Performance Profiler Report
===========================================
Process: 0 / 4
-------------------------------------------

Overall Timing:
  Total Time:       1.234 s
  Communication:    456 ms (36.98%)
  Computation:      778 ms (63.02%)
  Comm/Comp Ratio:  0.586
  Load Balance:     92.5%

Communication Statistics:
  Total Messages:   24
  Sends:           12
  Receives:        12
  Total Bytes:     96 KB
  Bytes Sent:      48 KB
  Bytes Received:  48 KB
  Avg Send Size:   4.00 KB
  Avg Recv Size:   4.00 KB
```

### CSV Format

Machine-readable format for plotting and analysis:

```csv
# MPI Performance Profiler CSV Report
# Process: 0 / 4
# Total Time: 1.234

# Timing Regions
Region,Calls,Total_Time,Avg_Time_Per_Call
"compute_step",100,0.778,0.00778
"exchange_ghost",4,0.456,0.114

# Communication by Tag
Tag,Sends,Receives,Bytes_Sent,Bytes_Received
"ghost_exchange",12,12.49152,49152
```

### JSON Format

Structured format for programmatic analysis:

```json
{
  "process_rank": 0,
  "num_processes": 4,
  "total_time": 1.234,
  "timing": {
    "communication_time": 0.456,
    "computation_time": 0.778,
    "idle_time": 0.0,
    "comm_percentage": 36.98,
    "comp_percentage": 63.02,
    "comm_comp_ratio": 0.586
  },
  "communication": {
    "total_messages": 24,
    "send_count": 12,
    "recv_count": 12,
    "total_bytes": 98304,
    "bytes_sent": 49152,
    "bytes_received": 49152
  }
}
```

## Visualization

### Using the Python Script

The profiling script automatically generates a Python visualization script:

```bash
# Run profiling
./scripts/profile_mpi.sh -p bin/benchmark_mpi_communication -n 4 -o results/

# Generate plots
python3 results/visualize.py results/
```

This will generate PNG plots showing:
- Communication time vs message size
- Throughput vs message size
- Other performance metrics

### Using Gnuplot

If gnuplot is available, a gnuplot script is also generated:

```bash
# Edit the script to customize plots
vim results/plot_timings.gnuplot

# Generate plots
gnuplot results/plot_timings.gnuplot
```

## Profiling Tools

### Built-in Profiling

Uses MPI_Wtime() for lightweight timing with minimal overhead.

### Score-P

A powerful instrumentation and measurement framework for parallel applications.

**Installation**: `sudo apt-get install scorep`

**Usage**:
```bash
./scripts/profile_mpi.sh -p bin/heat_equation -n 4 -t scorep -o scorep_results/
```

**Viewing traces**:
```bash
scorep-score --io=pretty scorep_results/scorep_results/profile.cubex
```

### Intel VTune Amplifier

Intel's advanced performance analysis tool.

**Installation**: Requires Intel oneAPI

**Usage**:
```bash
./scripts/profile_mpi.sh -p bin/heat_equation -n 4 -t vtune -o vtune_results/
```

**Viewing results**:
```bash
vtune -result vtune_results/vtune_results/
```

## Performance Metrics

### Load Balance Efficiency

Me how evenly work is distributed across processes:

```
Load Balance = Average Time / Maximum Time
```

A value of 1.0 indicates perfect load balance. Lower values indicate imbalance.

### Communication/Computation Ratio

```
C/C Ratio = Communication Time / Computation Time
```

Higher values indicate the application is communication-bound rather than compute-bound.

### Throughput

```
Throughput = Data Size / Time
```

Measured in MB/s for large messages.

### Latency

Measured in nanoseconds for small messages (ping-pong style).

## Best Practices

1. **Disable profiling in production builds**:
   ```cpp
   #ifdef ENABLE_PROFILING
       Profiler profiler(true, rank, size);
   #else
       Profiler profiler(false, rank, size);
   #endif
   ```

2. **Use warmup runs** to avoid measuring initialization overhead

3. **Profile realistic workloads**, not micro-benchmarks alone

4. **Check load balance** across all processes, not just the root

5. **Compare different MPI implementations** for your hardware

6. **Profile at multiple process counts** to identify scaling issues

## Troubleshooting

### Profiler shows all zeros

- Ensure profiler is enabled: `profiler.set_enabled(true);`
- Check that start/end pairs match
- Verify timing regions are properly nested

### Script fails to find MPI

- Ensure MPI is installed and in your PATH
- Check that `mpirun` command is available
- Verify MPI environment variables are set

### Score-P not instrumenting

- Ensure Score-P is installed
- Check that your MPI implementation is supported
- Verify the program is compiled with debugging symbols (-g)

### VTune results empty

- Ensure Intel VTune Amplifier is installed
- Check that you have permissions to run VTune
- Verify the Intel oneAPI environment is sourced

## Additional Resources

- MPI Standard: https://www.mpi-forum.org/docs/
- Score-P Documentation: https://scorep.github.io/
- Intel VTune: https://www.intel.com/content/www/us/en/developer/tools/oneapi/vtune-profiler.html
- MPI Performance Guide: https://www.mcs.anl.gov/research/projects/mpi/tutorial/

## File Reference

### Source Files

- `src/mpi/profiler.hpp` - Profiler class header
- `src/mpi/profiler.cpp` - Profiler class implementation

### Benchmark Files

- `tests/performance/benchmark_mpi_communication.cpp` - Comprehensive MPI communication benchmarks
- `tests/performance/benchmark_communication.cpp` - Communication pattern benchmarks
- `tests/performance/benchmark_overlap.cpp` - Communication/computation overlap benchmarks
- `tests/performance/benchmark_scaling.cpp` - Scaling benchmarks

### Script Files

- `scripts/profile_mpi.sh` - Automated Profiling script

### Example Files

- `examples/profiler_example.cpp` - Profiler usage example

## Summary

The MPI performance profiling tools provide:

1. **Easy-to-use API** for instrumenting your code
2. **Comprehensive benchmarks** for evaluating MPI performance
3. **Automated profiling** with external tools
4. **Multiple report formats** for analysis and visualization
5. **Minimal overhead** when profiling is disabled

For questions or issues, please refer to the project documentation or open an issue on GitHub.
