#!/bin/bash

# profile_mpi.sh - MPI Performance Profiling Script
#
# This script runs MPI programs with various profiling tools and collects
# performance metrics. It supports:
# - Native MPI profiling (MPI_Pcontrol)
# - Score-P instrumentation
# - Intel VTune Amplifier
# - Built-in timer profiling
#
# Usage: ./profile_mpi.sh -n 4 -o results/ -v -t none

set -e

# Default values
NUM_PROCS=4
OUTPUT_DIR="results"
VERBOSE=false
TOOL="none"
PROGRAM=""
PROGRAM_ARGS=""
WARMUP_RUNS=1
MEASURED_RUNS=3

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Print usage information
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

MPI Performance Profiling Script

OPTIONS:
    -n NUM_PROCES     Number of MPI processes (default: 4)
    -o OUTPUT_DIR     Output directory for results (default: results/)
    -p PROGRAM        MPI program to profile (required)
    -a "ARGS"         Arguments to pass to the program
    -t TOOL           Profiling tool: none, scorep, vtune (default: none)
    -w WARMUP         Number of warmup runs (default: 1)
    -m MEASURED       Number of measured runs (default: 3)
    -v                Verbose output
    -h                Show this help message

EXAMPLES:
    # Run with built-in profiling
    $0 -p ./heat_equation -n 4 -o my_results -v

    # Run with Score-P
    $0 -p ./heat_equation -n 4 -t scorep -o scorep_results

    # Run with Intel VTune
    $0 -p ./heat_equation -n 4 -tvtune -o vtune_results

EOF
    exit 1
}

# Parse command line arguments
while getopts "n:o:p:a:t:w:m:vh" opt; do
    case $opt in
        n) NUM_PROCS="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        p) PROGRAM="$OPTARG" ;;
        a) PROGRAM_ARGS="$OPTARG" ;;
        t) TOOL="$OPTARG" ;;
        w) WARMUP_RUNS="$OPTARG" ;;
        m) MEASURED_RUNS="$OPTARG" ;;
        v) VERBOSE=true ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check if program is specified
if [ -z "$PROGRAM" ]; then
    echo -e "${RED}Error: Program to profile not specified (-p).${NC}"
    usage
fi

# Check if program exists
if [ ! -f "$PROGRAM" ]; then
    echo -e "${RED}Error: Program '$PROGRAM' not found.${NC}"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Print configuration
echo -e "${BLUE}===========================================${NC}"
echo -e "${BLUE}MPI Performance Profiling${NC}"
echo -e "${BLUE}===========================================${NC}"
echo "Program: $PROGRAM"
echo "Arguments: $PROGRAM_ARGS"
echo "Processes: $NUM_PROCS"
echo "Tool: $TOOL"
echo "Output: $OUTPUT_DIR"
echo "Warmup runs: $WARMUP_RUNS"
echo "Measured runs: $MEASURED_RUNS"
echo -e "${BLUE}===========================================${NC}"
echo ""

# Check if MPI is available
if ! command -v mpirun &> /dev/null; then
    echo -e "${RED}Error: mpirun not found. Please ensure MPI is installed.${NC}"
    exit 1
fi

# Function to run MPI program
run_mpi() {
    local output_file="$1"
    shift

    if [ "$VERBOSE" = true ]; then
        echo -e "${GREEN}Running: mpirun -np $NUM_PROCES $PROGRAM $PROGRAM_ARGS${NC}"
    fi

    mpirun -np "$NUM_PROCS" "$PROGRAM" $PROGRAM_ARGS 2>&1 | tee "$output_file"
}

# Function to run with Score-P
run_scorep() {
    local output_file="$1"
    local scorep_dir="$OUTPUT_DIR/scorep_results"

    if [ "$VERBOSE" = true ]; then
        echo -e "${GREEN}Running with Score-P instrumentation${NC}"
    fi

    # Create Score-P directory
    mkdir -p "$scorep_dir"

    # Run with Score-P instrumentation
    scorep --nocompiler --user --output="$scorep_dir/profile" \
        mpirun -np "$NUM_PROCS" "$PROGRAM" $PROGRAM_ARGS 2>&1 | tee "$output_file"

    # Generate Score-P reports
    if [ -f "$scorep_dir/profile.cubex" ]; then
        echo -e "${GREEN}Score-P trace generated: $scorep_dir/profile.cubex${NC}"
        echo -e "${YELLOW}View with: scorep-score --io=pretty $scorep_dir/profile.cubex${NC}"
    fi
}

# Function to run with Intel VTune
run_vtune() {
    local output_file="$1"
    local vtune_dir="$OUTPUT_DIR/vtune_results"

    if [ "$VERBOSE" = true ]; then
        echo -e "${GREEN}Running with Intel VTune Amplifier${NC}"
    fi

    # Create VTune directory
    mkdir -p "$vtune_dir"

    # Check if VTune is available
    if ! command -v amplxe-cl &> /dev/null; then
        echo -e "${YELLOW}Warning: VTune Amplifier not found. Falling back to MPI profiling.${NC}"
        run_mpi "$output_file"
        return
    fi

    # Run with VTune
    amplxe-cl -collect hotspots -result-dir "$vtune_dir" -- \
        mpirun -np "$NUM_PROCS" "$PROGRAM" $PROGRAM_ARGS 2>&1 | tee "$output_file"

    echo -e "${GREEN}VTune results saved to: $vtune_dir${NC}"
}

# Main profiling loop
echo -e "${BLUE}Starting profiling...${NC}"
echo ""

# Warmup runs
echo -e "${YELLOW}Warmup runs...${NC}"
for ((i=1; i<=WARMUP_RUNS; i++)); do
    if [ "$VERBOSE" = true ]; then
        echo "Warmup run $i/$WARMUP_RUNS"
    fi

    case "$TOOL" in
        scorep)
            # Don't instrument warmup runs for Score-P
            run_mpi "$OUTPUT_DIR/warmup_$i.log"
            ;;
        vtune)
            # Don't profile warmup runs for VTune
            run_mpi "$OUTPUT_DIR/warmup_$i.log"
            ;;
        *)
            run_mpi "$OUTPUT_DIR/warmup_$i.log"
            ;;
    esac
done

echo ""
echo -e "${YELLOW}Measured runs...${NC}"

# Measured runs
for ((i=1; i<=MEASURED_RUNS; i++)); do
    echo -e "${GREEN}Run $i/$MEASURED_RUNS${NC}"

    local log_file="$OUTPUT_DIR/run_$i.log"

    case "$TOOL" in
        scorep)
            run_scorep "$log_file"
            ;;
        vtune)
            run_vtune "$log_file"
            ;;
        none)
            run_mpi "$log_file"
            ;;
        *)
            echo -e "${YELLOW}Unknown tool: $TOOL. Using MPI profiling.${NC}"
            run_mpi "$log_file"
            ;;
    esac

    echo ""
done

# Generate summary
echo -e "${BLUE}===========================================${NC}"
echo -e "${BLUE}Profiling Complete${NC}"
echo -e "${BLUE}===========================================${NC}"

# Extract timing information if available
echo -e "\n${YELLOW}Timing Summary:${NC}"
for ((i=1; i<=MEASURED_RUNS; i++)); do
    if [ -f "$OUTPUT_DIR/run_$i.log" ]; then
        echo -n "Run $i: "
        grep -i "time" "$OUTPUT_DIR/run_$i.log" | tail -1 || echo "No timing info found"
    fi
done

# Generate simple visualization script if gnuplot is available
if command -v gnuplot &> /dev/null && [ "$TOOL" = "none" ]; then
    echo -e "\n${YELLOW}Generating visualization...${NC}"

    # Create gnuplot script
    cat > "$OUTPUT_DIR/plot_timings.gnuplot" << 'EOF'
set terminal png size 800,600
set output 'timings_plot.png'
set title 'MPI Performance Timings'
set xlabel 'Run Number'
set ylabel 'Time (s)'
set grid

# Parse timing data from logs (example - adjust based on actual output format)
# This is a placeholder - actual parsing depends on program output
plot '-' w lp title 'Execution Time'

# Add data points here after parsing logs
EOF

    echo -e "${YELLOW}Gnuplot script generated: $OUTPUT_DIR/plot_timings.gnuplot${NC}"
    echo -e "${YELLOW}Edit the script and run: gnuplot $OUTPUT_DIR/plot_timings.gnuplot${NC}"
fi

# Alternative: Generate Python visualization script
cat > "$OUTPUT_DIR/visualize.py" << 'EOF'
#!/usr/bin/env python3
"""
Simple Python script to visualize MPI profiling results.
Requires: pandas, matplotlib
"""

import os
import sys
import glob

try:
    import pandas as pd
    import matplotlib.pyplot as plt
except ImportError:
    print("Error: pandas and matplotlib are required.")
    print("Install with: pip install pandas matplotlib")
    sys.exit(1)

def visualize_results(output_dir):
    """Generate visualization from profiling results."""

    # Look for CSV files
    csv_files = glob.glob(os.path.join(output_dir, "*.csv"))

    if not csv_files:
        print("No CSV files found in output directory.")
        return

    # Process each CSV file
    for csv_file in csv_files:
        try:
            df = pd.read_csv(csv_file)

            # Skip comment lines if present
            df = df[~df.iloc[:, 0].str.startswith('#', na=False)]

            # Create plot
            plt.figure(figsize=(10, 6))

            if 'avg_time_s' in df.columns and 'message_size' in df.columns:
                plt.plot(df['message_size'], df['avg_time_s'], marker='o')
                plt.xlabel('Message Size')
                plt.ylabel('Average Time (s)')
                plt.title('Communication Time vs Message Size')
                plt.xscale('log')
                plt.yscale('log')
                plt.grid(True)
            elif 'throughput_mb_s' in df.columns:
                plt.plot(df['message_size'], df['throughput_mb_s'], marker='o')
                plt.xlabel('Message Size')
                plt.ylabel('Throughput (MB/s)')
                plt.title('Throughput vs Message Size')
                plt.xscale('log')
                plt.grid(True)
            else:
                # Generic plot
                df.plot(kind='bar')
                plt.title('Performance Metrics')
                plt.xticks(rotation=45)

            # Save plot
            basename = os.path.splitext(os.path.basename(csv_file))[0]
            plot_file = os.path.join(output_dir, f"{basename}.png")
            plt.tight_layout()
            plt.savefig(plot_file)
            plt.close()

            print(f"Plot saved: {plot_file}")

        except Exception as e:
            print(f"Error processing {csv_file}: {e}")

    print("\nVisualization complete!")
    print(f"Plots saved to: {output_dir}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        output_dir = sys.argv[1]
    else:
        output_dir = "."

    visualize_results(output_dir)
EOF

chmod +x "$OUTPUT_DIR/visualize.py"
echo -e "${GREEN}Python visualization script: $OUTPUT_DIR/visualize.py${NC}"
echo -e "${YELLOW}Run with: python3 $OUTPUT_DIR/visualize.py $OUTPUT_DIR${NC}"

# Print final instructions
echo ""
echo -e "${BLUE}===========================================${NC}"
echo -e "${BLUE}Results Summary${NC}"
echo -e "${BLUE}===========================================${NC}"
echo "Output directory: $OUTPUT_DIR"
echo ""
echo "Files available:"
ls -lh "$OUTPUT_DIR"
echo ""
echo "Next steps:"
echo "1. Review the log files for timing information"
echo "2. Run: python3 $OUTPUT_DIR/visualize.py $OUTPUT_DIR"
echo "3. (Optional) Edit plot_timings.gnuplot and run gnuplot for custom plots"
echo ""
