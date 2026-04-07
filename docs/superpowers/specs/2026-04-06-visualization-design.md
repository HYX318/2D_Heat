# 2D Heat Equation Solver - Python Visualization Design

**Date:** 2026-04-06
**Purpose:** Add end-to-end Python visualization that runs the MPI solver and automatically visualizes plots/animations

---

## Overview

Create a Python visualization system that provides a complete workflow: run the compiled MPI solver, monitor output file generation, and automatically visualize the results as heat maps or animations.

---

## Requirements

### Functional Requirements

1. **Automatic Solver Execution**
   - Run compiled `mpirun -n X ./bin/heat_equation` from Python
   - Support configurable solver parameters (grid size, time steps, solver type, etc.)
   - Handle both serial and MPI execution

2. **Preset Configuration**
   - Provide predefined configurations (small, medium, large)
   - Allow command-line overrides of preset values
   - YAML-based configuration file

3. **Output File Detection**
   - Automatically detect output file format (Text or VTK)
   - Monitor file generation during solver execution
   - Collect all generated output files

4. **Visualization**
   - Plot single time step as 2D heat map
   - Create animation from multiple time steps
   - Auto-detect whether to show plot or animation based on file count
   - Display matplotlib windows
   - Save images to files

5. **CLI Interface**
   - Simple command-line interface
   - Support preset selection
   - Support parameter overrides
   - Clear help documentation

### Non-Functional Requirements

- **Dependencies**: matplotlib, numpy, pyyaml, tqdm
- **Python Version**: 3.7+
- **Error Handling**: Graceful handling of missing executables, file errors
- **Performance**: Non-blocking file monitoring with progress indication

---

## Architecture

```
visualization/
├── run_and_visualize.py    # Main entry point - CLI and orchestration
├── requirements.txt          # Python dependencies
├── config/
│   └── presets.yaml         # Preset solver configurations
├── core/
│   ├── __init__.py
│   ├── solver_runner.py     # MPI solver execution
│   ├── solution_reader.py   # File parsing (Text/VTK)
│   └── output_monitor.py    # File monitoring utilities
└── renderers/
    ├── __init__.py
    ├── heatmap.py           # Single frame rendering
    └── animation.py        # Multi-frame animation
```

---

## Components

### 1. `run_and_visualize.py` (Main Entry Point)

**Responsibilities:**
- Parse command-line arguments
- Load preset configurations
- Merge solver parameters
- Orchestrate the workflow: run → monitor → visualize

**Key Functions:**
```python
def main():
    """Main entry point with CLI argument parsing"""
    args = parse_arguments()
    config = load_and_merge_config(args)
    run_workflow(config)

def run_workflow(config: Dict):
    """Execute complete workflow: run solver, monitor, visualize"""
    runner = SolverRunner(config)
    files = runner.run_and_collect()
    visualize_results(files, config)
```

---

### 2. `core/solver_runner.py`

**Class:** `SolverRunner`

**Responsibilities:**
- Construct mpirun command with all parameters
- Execute solver as subprocess
- Monitor output file generation
- Return list of generated files

**Key Methods:**
```python
class SolverRunner:
    def __init__(self, config: Dict, build_dir: str = "../build"):
        """Initialize with solver configuration and build directory"""
    
    def build_command(self) -> List[str]:
        """Construct full mpirun command"""
        # Example: ["mpirun", "-n", "4", "./bin/heat_equation", "-nx", "100", ...]
    
    def run(self) -> subprocess.Popen:
        """Execute solver and return process handle"""
    
    def wait_for_completion(self, process: subprocess.Popen) -> int:
        """Wait for solver to finish and monitor progress"""
    
    def collect_output_files(self) -> List[str]:
        """Find all generated output files matching prefix"""
    
    def run_and_collect(self) -> List[str]:
        """Run solver and collect output files in one call"""
```

---

### 3. `core/solution_reader.py`

**Class:** `SolutionReader`

**Responsibilities:**
- Detect file format (Text or VTK)
- Parse both formats into standardized format
- Extract metadata (grid size, time, step)

**Key Methods:**
```python
class SolutionReader:
    @staticmethod
    def detect_format(filename: str) -> str:
        """Detect if file is Text or VTK format"""
    
    @staticmethod
    def read_text_file(filename: str) -> Dict:
        """
        Parse Text format output
        Returns: {
            'x': ndarray, 'y': ndarray, 'temperature': ndarray,
            'nx': int, 'ny': int, 'time': float, 'step': int
        }
        """
    
    @staticmethod
    def read_vtk_file(filename: str) -> Dict:
        """Parse VTK format and return same structure as read_text_file"""
    
    @staticmethod
    def read_file(filename: str) -> Dict:
        """Auto-detect format and read file"""
    
    @staticmethod
    def read_pattern(pattern: str) -> List[Dict]:
        """Read all files matching glob pattern, return sorted by step"""
```

---

### 4. `core/output_monitor.py`

**Functions:**
```python
def monitor_directory(directory: str, pattern: str, interval: float = 0.5) -> Iterator[str]:
    """
    Generator that yields new files as they appear
    Useful for showing progress during long solver runs
    """

def wait_for_files(expected_count: int, directory: str, pattern: str, timeout: int = 300) -> List[str]:
    """
    Wait until expected number of files exist
    Returns list of file paths when ready
    """
```

---

### 5. `renderers/heatmap.py`

**Functions:**
```python
def plot_heat_map(
    x: ndarray,
    y: ndarray,
    temperature: ndarray,
    title: str = "Heat Distribution",
    colormap: str = "hot",
    show: bool = True,
    save_path: str = None,
    colorbar: bool = True
) -> Figure:
    """
    Plot 2D heat distribution as contourf or pcolormesh
    Shows grid colorbar by default
    Returns matplotlib Figure for further customization
    """
```

---

### 6. `renderers/animation.py`

**Functions:**
```python
def create_animation(
    files_data: List[Dict],  # List of parsed solution data
    output_path: str = None,
    show: bool = True,
    fps: int = 10,
    colormap: str = "hot"
) -> FuncAnimation:
    """
    Create animation from multiple time steps
    Returns matplotlib FuncAnimation object
    """
```

---

## Configuration

### `config/presets.yaml`

```yaml
presets:
  small:
    processes: 2
    nx: 50
    ny: 50
    nt: 100
    dt: 0.01
    solver: Jacobi
    scheme: ImplicitEuler
    output_prefix: "output/small"
    
  medium:
    processes: 4
    nx: 100
    ny: 100
    nt: 500
    dt: 0.005
    solver: SOR
    scheme: CrankNicolson
    output_prefix: "output/medium"
    
  large:
    processes: 8
    nx: 200
    ny: 200
    nt: 2000
    dt: 0.002
    solver: CG
    scheme: CrankNicolson
    output_prefix: "output/large"
```

### Default Parameter Mapping

| YAML Key | Solver CLI Arg | Type | Default |
|-----------|----------------|-------|----------|
| processes | mpirun -n | int | 1 |
| nx | -nx | int | 100 |
| ny | -ny | int | 100 |
| nt | -nt | int | 100 |
| dt | -dt | float | 0.01 |
| solver | -solver | str | Jacobi |
| scheme | -scheme | str | ImplicitEuler |
| output_prefix | -o | str | output/solution |

---

## CLI Interface

### Command-Line Arguments

```bash
python run_and_visualize.py [OPTIONS]

Preset Selection:
  --preset PRESET       Use preset config (small/medium/large)
  
Solver Parameters:
  -n, --processes N    Number of MPI processes
  -nx N                Grid points in x-direction
  -ny N                Grid points in y-direction
  -nt N                Number of time steps
  -dt FLOAT            Time step size
  --solver TYPE         Solver type (Jacobi/SOR/CG)
  --scheme TYPE         Time scheme (ImplicitEuler/CrankNicolson)
  
Output Options:
  - -o, --output PREFIX   Output file prefix
  --format FMT        Output format (Text/VTK)
  --build-dir PATH     Build directory (default: ../build)
  
Visualization (automatic):
  --no-plot           Skip plotting
  --save-only         Save images, don't show
  --dpi N             Image resolution for saved files
  
Other:
  -h, --help          Show help message
  -v, --verbose       Verbose output
```

### Usage Examples

```bash
# Use preset
python run_and_visualize.py --preset medium

# Preset with overrides
python run_and_visualize.py --preset small --nt 500 --solver CG

# Full custom run
python run_and_visualize.py \
    -n 4 -nx 200 -ny 200 -nt 1000 \
    -solver SOR -scheme CrankNicolson \
    -output results/my_run

# Custom build directory
python run_and_visualize.py --preset large --build-dir /path/to/build
```

---

## Workflow

### Complete Execution Flow

```
1. Parse command-line arguments
   ↓
2. Load preset configuration (if --preset specified)
   ↓
3. Merge command-line overrides (CLI args take precedence)
   ↓
4. Validate configuration
   ↓
5. Construct mpirun command with all parameters
   ↓
6. Execute solver as subprocess
   ↓
7. Monitor output directory for new files (with progress indicator)
   ↓
8. Wait for solver completion
   ↓
9. Collect all generated output files
   ↓
10. Sort files by time step
    ↓
11. If single file:
        - Parse file
        - Plot as heat map
        - Display window
    If multiple files:
        - Parse all files
        - Create animation
        - Display animation
    ↓
12. Save images/animations if requested
    ↓
13. Exit
```

---

## Error Handling

### Scenarios and Responses

| Scenario | Response |
|----------|----------|
| Build directory not found | Error: "Build directory not found. Run 'mkdir build && cd build && cmake .. && make' first" |
| heat_equation executable not found | Error: "Solver executable not found. Ensure build is complete" |
| No output files generated | Warning: "No output files found. Check solver ran successfully" |
| Unknown file format | Error: "Unsupported output format: .xxx" |
| Preset not found | Error: "Preset 'xxx' not found. Available: small, medium, large" |
| Invalid parameter | Error: "Invalid value for --nt: must be positive integer" |
| Solver execution failed | Error: "Solver exited with code N. Check stderr for details" |

---

## Testing Strategy

### Unit Tests

1. **Solution Reader Tests**
   - Test Text format parsing
   - Test VTK format parsing
   - Test format detection

2. **Command Builder Tests**
   - Verify mpirun command construction
   - Test parameter mapping

3. **Configuration Tests**
   - Test preset loading
   - Test parameter merging
   - Test default values

### Integration Tests

1. **End-to-End Workflow**
   - Run solver with preset
   - Verify output files generated
   - Verify plots created

2. **Visualization Tests**
   - Test heatmap rendering
   - Test animation creation (with reduced frames)

---

## Future Enhancements

1. **3D Visualization** - Add 3D surface plots
2. **Real-time Monitoring** - Live updating plot during solver execution
3. **Comparison Mode** - Overlay multiple solver results
4. **Export Formats** - Support exporting as GIF, MP4
5. **Interactive Plotting** - Add hover tooltips, zoom controls
6. **Performance Metrics** - Plot convergence history, timing data

---

## Dependencies

**requirements.txt:**
```
matplotlib>=3.5.0
numpy>=1.21.0
pyyaml>=6.0
tqdm>=4.62.0
```

---

## Notes

- Script assumes it's run from `visualization/` directory
- Relative paths used for build directory access
- MPI must be installed and available in PATH
- Output files are expected to match prefix + pattern
