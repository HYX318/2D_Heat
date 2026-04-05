# MPI Communication Performance Benchmarks

This directory contains comprehensive benchmarks for evaluating MPI communication performance in the 2D Heat Equation solver.

## Overview

The benchmark suite includes three main categories of tests:

1. **Communication Pattern Benchmarks** (`benchmark_communication.cpp`)
   - Ping-Pong latency and bandwidth tests
   - Allreduce collective operation benchmarks
   - Broadcast efficiency tests
   - Barrier synchronization overhead

2. **Scaling Tests** (`benchmark_scaling.cpp`)
   - Strong scaling: Fixed total workload, increasing processors
   - Weak scaling: Fixed per-processor workload, increasing processors

3. **Communication/Computation Overlap** (`benchmark_overlap.cpp`)
   - Non-blocking communication with compute overlap
   - Comparison of async vs sync performance
   - Overlap ratio analysis

## File Structure

```
tests/performance/
├── benchmark_communication.cpp     # Communication pattern benchmarks
├── benchmark_scaling.cpp            # Strong/weak scaling tests
├── benchmark_overlap.cpp            # Communication/computation overlap tests
├── CMakeLists.txt                   # Build configuration
└── README.md                        # This file
```

## Building the Benchmarks

### Prerequisites

- CMake 3.16 or higher
- C++17 compatible compiler
- MPI implementation (OpenMPI, MPICH, etc.)

### Build Instructions

```bash
# Navigate to build directory
cd /path/to/2D_Heat/build

# Configure CMake
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build all benchmarks
make benchmark_communication benchmark_scaling benchmark_overlap

# Or build all targets
make
```

### Build Output

Executables are placed in `build/bin/`:
- `benchmark_communication`
- `benchmark_scaling`
- `benchmark_overlap`

## Running the Benchmarks

### Communication Benchmarks

```bash
# Run all communication benchmarks with 4 processes
mpirun -np 4 ./bin/benchmark_communication --all

# Run specific benchmark
mpirun -np 4 ./bin/benchmark_communication --test pingpong
mpirun -np 4 ./bin/benchmark_communication --test allreduce
mpirun -np 4 ./bin/benchmark_communication --test broadcast
mpirun -np 4 ./bin/benchmark_communication --test barrier

# Customize iterations and output
mpirun -np 4 ./bin/benchmark_communication --all --iterations 5000 --output my_results.csv
```

#### Communication Benchmark Options

- `--all`: Run all benchmarks (default)
- `--test <name>`: Run specific test (pingpong, allreduce, broadcast, barrier)
- `--iterations <n>`: Number of iterations per test (default: 1000)
- `--output <file>`: CSV output filename (default: benchmark_results.csv)
- `--help`: Show help message

### Scaling Benchmarks

```bash
# Run strong scaling test
mpirun -np 4 ./bin/benchmark_scaling --type strong --strong-size 1024

# Run weak scaling test
mpirun -np 4.0 ./bin/bb/benchmark_scaling --type weak --weak-per-proc 256

# Run both strong and weak scaling
mpirun -np 4 ./bin/benchmark_scaling --type both

# Generate scaling curve (run multiple times)
for np in 1 2 4 8 16; do
    mpirun -np $np ./bin/benchmark_scaling --type both --output scaling_np${np}.csv
done
```

#### Scaling Benchmark Options

- `--type <type>`: Test type (strong, weak, both) (default: both)
- `--max <n>`: Maximum number of processes (default: current)
- `--strong-size <n>`: Total grid size for strong scaling (default: 1024)
- `--weak-per-proc <n>`: Per-process grid size for weak scaling (default: 256)
- `--compute-iter <n>`: Number of compute iterations (default: 100)
- `--comm-iter <n>`: Number of communication iterations (default: 1000)
- `--output <file>`: CSV output filename (default: scaling_results.csv)
- `--help`: Show help message

### Overlap Benchmarks

```bash
# Compare async vs sync performance
mpirun -np 4 ./bin/benchmark_overlap --test compare

# Run only async overlap benchmark
mpirun -np 4 ./bin/benchmark_overlap --test async

# Run only sync benchmark
mpirun -np 4 ./bin/benchmark_overlap --test sync

# Customize grid size and message size
mpirun -np 4 ./bin/benchmark_overlap --test compare --grid-size 512 --message-size 2048
```

#### Overlap Benchmark Options

- `--test <type>`: Test type (async, sync, compare) (default: compare)
- `--grid-size <n>`: Grid size (default: 256)
- `--message-size <n>`: Message size in doubles (default: 1024)
- `--iterations <n>`: Number of iterations (default: 100)
- `--output <file>`: CSV output filename (default: overlap_results.csv)
- `--help`: Show help message

## Output Format

### Console Output

Benchmarks produce formatted console output:

```
========================================
   MPI Communication Performance Benchmark
========================================
Number of processes: 4
Iterations per test: 1000
CSV output: benchmark_results.csv

=== Ping-Pong Benchmark (Send/Recv) ===
Iterations: 1000

1.0B    :        12.34 µs (latency)     162.34 MB/s (bandwidth)
1.0KB   :        45.67 µs (latency)      44.32 MB/s (bandwidth)
1.0MB   :      1234.56 µs (latency)       1.64 MB/s (bandwidth)
10.0MB  :     12345.67 µs (latency)       1.64 MB/s (bandwidth)
```

### CSV Output

Benchmarks also produce CSV files for analysis:

**Communication Benchmark CSV**:
```csv
test_type,processes,message_size,time_us,bandwidth_mbps,operation
pingpong_sendrecv,2,1,12.34,162.34,round-trip
pingpong_sendrecv,2,1024,45.67,44.32,round-trip
allreduce,4,2048,123.45,33.12,SUM
```

**Scaling Benchmark CSV**:
```csv
test_type,processes,problem_size,per_proc_size,time_ms,efficiency,gflops
strong_scaling,4,1048576,262144,123.45,95.67,3.45
weak_scaling,4,65536,16384,45.67,98.23,1.23
```

**Overlap Benchmark CSV**:
```csv
test_type,processes,message_size,iterations,time_ms,compute_time_ms,comm_time_ms,overlap_ratio
async_overlap,4,1024,100,123.45,100.00,50.00,1.22
sync_no_overlap,4,1024,100,150.00,100.00,50.00,1.00
```

## Performance Metrics

### Latency

- Measured in microseconds (µs)
- Time for a single message to travel between processes
- Includes MPI overhead and network latency

### Bandwidth

- Measured in MB/s or GB/s
- Rate of data transfer between processes
- Calculated as: `message_size / time`

### Efficiency

- Percentage of ideal speedup
- Strong scaling: `(T1 / Tp) / p * 100`
- Weak scaling: Should remain near 100% as processors increase

### Overlap Ratio

- Metric for communication/computation overlap
- Values > 1.0 indicate good overlap
- Calculated as: `(compute_time + comm) / total_time`

### GFLOPS

- Computational performance
- Measured in billions of floating point operations per second
- Useful for comparing compute-bound scenarios

## Performance Analysis

### Strong Scaling Analysis

Strong scaling keeps the total problem size fixed while increasing the number of processors:

```bash
# Run strong scaling with different processor counts
for np in 1 2 4 8 16; do
    mpirun -np $np ./bin/benchmark_scaling \
        --type strong \
        --strong-size 1024 \
        --output strong_np${np}.csv
done

# Analyze results
# - Speedup should increase with more processors
# - Efficiency may decrease due to communication overhead
```

### Weak Scaling Analysis

Weak scaling keeps the per-processor workload fixed while increasing the number of processors:

```bash
# Run weak scaling with different processor counts
for np in 1 2 4 8 16; do
    mpirun -np $np ./bin/benchmark_scaling \
        --type weak \
        --weak-per-proc 256 \
        --output weak_np${np}.csv
done

# Analyze results
# - Total time should remain constant
# - Efficiency should stay near 100%
```

### Communication Overlap Analysis

Compare async vs sync performance to measure overlap effectiveness:

```bash
# Run comparison
mpirun -np 4 ./bin/benchmark_overlap --test compare

# Analyze overlap_ratio:
# - Values > 1.0 indicate effective overlap
# - Higher values mean better overlap
# - Ideal: overlap_ratio = (compute_time + comm_time) / max(compute_time, comm_time)
```

## Troubleshooting

### MPI Not Found

If you encounter `MPI not found` errors:

```bash
# Set MPI environment variables
export MPI_CXX=<your-mpicxx-path>
export MPI_ROOT=<your-mpi-install-path>

# Reconfigure CMake
cmake .. -DCMAKE_BUILD_TYPE=Release
```

### Build Errors

If you encounter build errors:

1. Check that MPI is properly installed
2. Verify C++17 support
3. Ensure all dependencies are available
4. Check CMake configuration messages

### Runtime Errors

If you encounter runtime errors:

1. Ensure mpirun is in your PATH
2. Verify MPI launcher is configured correctly
3. Check that the number of processes doesn't exceed available resources
4. Ensure network settings allow MPI communication

## Advanced Usage

### Integration with Profiling Tools

The benchmarks can be used with profiling tools:

```bash
# Intel VTune
mpirun -np 4 vtune -collect hotspots ./bin/benchmark_communication --all

# Score-P
scorep mpirun -np 4 ./bin/benchmark_communication --all

# perf (Linux)
mpirun -np 4 perf stat ./bin/benchmark_communication --all
```

### Custom Benchmarking

You can modify the benchmark sources to:

- Add custom communication patterns
- Test specific message sizes
- Implement domain-specific workloads
- Add custom metrics and output formats

## Contributing

When adding new benchmarks:

1. Follow the existing naming conventions
2. Include CSV output support
3. Add command-line parameter handling
4. Document the benchmark in this README
5. Update CMakeLists.txt if needed

## References

- [MPI Standard](https://www.mpi-forum.org/docs/)
- [OpenMPI Performance Tuning](https://www.open-mpi.org/doc/)
- [MPICH Performance Guide](https://www.mpich.org/guide/)
- [MPI Tutorial](https://mpitutorial.com/)

## License

These benchmarks are part of the 2D Heat Equation solver project.
See the main project LICENSE file for details.
