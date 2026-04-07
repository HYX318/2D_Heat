# Python Visualization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create an end-to-end Python visualization system that runs the MPI solver and automatically plots/animations the results

**Architecture:** Modular Python package with separate components for solver execution, file parsing, and rendering. CLI-driven workflow with preset configuration support.

**Tech Stack:** Python 3.7+, matplotlib, numpy, pyyaml, tqdm, subprocess

---

## File Structure

```
visualization/
├── requirements.txt          # Python dependencies
├── config/
│   └── presets.yaml         # Preset solver configurations
├── core/
│   ├── __init__.py         # Package marker
│   ├── solution_reader.py   # File parsing (Text/VTK formats)
│   ├── output_monitor.py    # File monitoring utilities
│   └── solver_runner.py     # MPI solver execution
├── renderers/
│   ├── __init__.py         # Package marker
│   ├── heatmap.py           # Single frame rendering
│   └── animation.py        # Multi-frame animation
└── run_and_visualize.py    # Main entry point - CLI and orchestration
```

---

### Task 1: Create visualization directory structure and requirements.txt

**Files:**
- Create: `visualization/requirements.txt`

- [ ] **Step 1: Create visualization directory and requirements.txt**

```bash
mkdir -p visualization/config visualization/core visualization/renderers
```

Create `visualization/requirements.txt`:
```txt
matplotlib>=3.5.0
numpy>=1.21.0
pyyaml>=6.0
tqdm>=4.62.0
```

- [ ] **Step 2: Commit**

```bash
git add visualization/requirements.txt
git commit -m "feat: add visualization directory structure and dependencies"
```

---

### Task 2: Create preset configuration file

**Files:**
- Create: `visualization/config/presets.yaml`

- [ ] **Step 1: Create presets.yaml with small, medium, large configurations**

Create `visualization/config/presets.yaml`:
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

- [ ] **Step 2: Commit**

```bash
git add visualization/config/presets.yaml
git commit -m "feat: add preset configurations for solver"
```

---

### Task 3: Create core package __init__.py

**Files:**
- Create: `visualization/core/__init__.py`

- [ ] **Step 1: Create empty __init__.py file**

Create `visualization/core/__init__.py`:
```python
"""
Core utilities for visualization package.
Handles solver execution, file parsing, and monitoring.
"""
```

- [ ] **Step 2: Commit**

```bash
git add visualization/core/__init__.py
git commit -m "feat: add core package initialization"
```

---

### Task 4: Create solution_reader.py with Text format parsing

**Files:**
- Create: `visualization/core/solution_reader.py`
- Test: Manual testing with sample data

- [ ] **Step 1: Create solution_reader.py with SolutionReader class and Text format parser**

Create `visualization/core/solution_reader.py`:
```python
"""
Solution file reader supporting Text and VTK formats.
"""
import numpy as np
import os
from typing import Dict, List


class SolutionReader:
    """Reads solution files from heat equation solver."""

    @staticmethod
    def detect_format(filename: str) -> str:
        """Detect if file is Text or VTK format."""
        with open(filename, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('# vtk DataFile Version'):
                return 'vtk'
            else:
                return 'text'

    @staticmethod
    def read_text_file(filename: str) -> Dict:
        """
        Parse Text format output from solver.
        Returns: {
            'x': ndarray, 'y': ndarray, 'temperature': ndarray,
            'nx': int, 'ny': int, 'time': float, 'step': int
        }
        """
        data = []
        nx = 0
        ny = 0
        time = 0.0
        step = 0
        lx = 1.0
        ly = 1.0

        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#'):
                    # Parse metadata
                    if 'time' in line:
                        # Extract time value from line like "# Solution at time = 1.5"
                        parts = line.split('=')
                        if len(parts) > 1:
                            time = float(parts[1].split(',')[0].strip())
                    elif 'Grid' in line:
                        # Extract nx and ny from line like "# Grid: 100 x 100"
                        parts = line.split(':')
[1].strip()
                        nx_str, ny_str = grid_info.split('x')
                        nx = int(nx_str.strip())
                        ny = int(ny_str.strip())
                    elif 'Physical domain' in line:
                        # Extract lx and ly
                        parts = line.split('[')[1].split(']')[0]
                        lx_str, ly_str = parts.split(',')
                        lx = float(lx_str.strip())
                        ly = float(ly_str.strip())
                    elif 'step' in line and 'time' not in line:
                        # Extract step from line like "# Solution at time = 0, step = 0"
                        parts = line.split('=')
                        if len(parts) > 2:
                            step = int(parts[2].strip())
                elif line and not line.startswith('#'):
                    # Parse data line: I J x y u(x,y)
                    parts = line.split()
                    if len(parts) >= 5:
                        i, j, x, y, u = int(parts[0]), int(parts[1]), \
                                              float(parts[2]), float(parts[3]), float(parts[4])
                        data.append((i, j, x, y, u))

        # Convert to numpy arrays
        data = np.array(data)
        x = data[:, 2].reshape(ny, nx)
        y = data[:, 3].reshape(ny, nx)
        temperature = data[:, 4].reshape(ny, nx)

        return {
            'x': x,
            'y': y,
            'temperature': temperature,
            'nx': nx,
            'ny': ny,
            'time': time,
            'step': step,
            'lx': lx,
            'ly': ly
        }

    @staticmethod
    def read_file(filename: str) -> Dict:
        """Auto-detect format and read file."""
        format_type = SolutionReader.detect_format(filename)
        if format_type == 'vtk':
            return SolutionReader.read_vtk_file(filename)
        else:
            return SolutionReader.read_text_file(filename)

    @staticmethod
    def read_vtk_file(filename: str) -> Dict:
        """
        Parse VTK format and return same structure as read_text_file.
        """
        points = []
        temperatures = []
        nx = 0
        ny = 0
        time = 0.0
        step = 0

        with open(filename, 'r') as f:
            in_points_section = False
            in_scalar_section = False
            reading_data = False

            for line in f:
                line = line.strip()

                # Extract time from header
                if 'time = ' in line and 'step = ' in line:
                    time_part = line.split('time = ')[1].split(',')[0]
                    step_part = line.split('step = ')[1]
                    time = float(time_part)
                    step = int(step_part)

                # Parse DIMENSIONS for nx, ny
                if line.startswith('DIMENSIONS'):
                    dims = line.split()[1:3]
                    nx = int(dims[0])
                    ny = int(dims[1])

                # Enter POINTS section
                if line == 'POINTS':
                    continue
                if in_points_section and not line.startswith('POINT_DATA'):
                    if not line.startswith('POINTS'):
                        point_data = line.split()
                        if len(point_data) >= 2:
                            points.append((float(point_data[0]), float(point_data[1])))

                # Start POINTS section
                if line.startswith('POINTS'):
                    in_points_section = True
                    continue

                # Exit POINTS section
                if line.startswith('POINT_DATA'):
                    in_points_section = False
                    continue

                # Enter SCALARS section
                if line.startswith('SCALARS'):
                    in_scalar_section = True
                    continue

                # Skip LOOKUP_TABLE line
                if in_scalar_section and line.startswith('LOOKUP_TABLE'):
                    reading_data = True
                    continue

                # Read temperature data
                if in_scalar_section and reading_data:
                    if line:
                        temperatures.append(float(line))

        # Convert to arrays
        points = np.array(points)
        x = points[:, 0].reshape(ny, nx)
        y = points[:, 1].reshape(ny, nx)
        temperature = np.array(temperatures).reshape(ny, nx)

        return {
            'x': x,
            'y': y,
            'temperature': temperature,
            'nx': nx,
            'ny': ny,
            'time': time,
            'step': step,
            'lx': x.max() if x.size > 0 else 1.0,
            'ly': y.max() if y.size > 0 else 1.0
        }

    @staticmethod
    def read_pattern(pattern: str) -> List[Dict]:
        """Read all files matching glob pattern, return sorted by step."""
        import glob

        files = sorted(glob.glob(pattern))
        results = []

        for filename in files:
            try:
                data = SolutionReader.read_file(filename)
                results.append(data)
            except Exception as e:
                print(f"Warning: Could not read {filename}: {e}")

        # Sort by step if available
        results.sort(key=lambda x: x.get('step', 0))

        return results
```

- [ ] **Step 2: Commit**

```bash
git add visualization/core/solution_reader.py
git commit -m "feat: add solution reader for Text and VTK formats"
```

---

### Task 5: Create output_monitor.py

**Files:**
- Create: `visualization/core/output_monitor.py`

- [ ] **Step 1: Create output_monitor.py with file monitoring utilities**

Create `visualization/core/output_monitor.py`:
```python
"""
Output file monitoring utilities.
"""
import os
import glob
import time
from typing import Iterator, List


def monitor_directory(directory: str, pattern: str, interval: float = 0.5) -> Iterator[str]:
    """
    Generator that yields new files as they appear.
    Useful for showing progress during long solver runs.
    """
    if not os.path.exists(directory):
        os.makedirs(directory)

    seen_files = set()

    while True:
        current_files = set(glob.glob(os.path.join(directory, pattern)))
        new_files = current_files - seen_files

        for filename in sorted(new_files):
            yield filename

        seen_files = current_files
        time.sleep(interval)


def wait_for_files(expected_count: int, directory: str, pattern: str, timeout: int = 300) -> List[str]:
    """
    Wait until expected number of files exist.
    Returns list of file paths when ready.
    """
    if not os.path.exists(directory):
        os.makedirs(directory)

    start_time = time.time()

    while time.time() - start_time < timeout:
        files = glob.glob(os.path.join(directory, pattern))
        if len(files) >= expected_count:
            return sorted(files)
        time.sleep(0.5)

    # Timeout reached
    current_files = glob.glob(os.path.join(directory, pattern))
    return sorted(current_files)
```

- [ ] **Step 2: Commit**

```bash
git add visualization/core/output_monitor.py
git commit -m "feat: add output file monitoring utilities"
```

---

### Task 6: Create solver_runner.py

**Files:**
- Create: `visualization/core/solver_runner.py`

- [ ] **Step 1: Create solver_runner.py with SolverRunner class**

Create `visualization/core/solver_runner.py`:
```python
"""
MPI solver execution and management.
"""
import subprocess
import os
import glob
from typing import Dict, List
from .output_monitor import monitor_directory


class SolverRunner:
    """Manages execution of the MPI heat equation solver."""

    # Parameter mapping from config keys to solver CLI arguments
    PARAM_MAPPING = {
        'nx': '-nx',
        'ny': '-ny',
        'nt': '-nt',
        'dt': '-dt',
        'solver': '-solver',
        'scheme': '-scheme',
        'output_prefix': '-o',
    }

    def __init__(self, config: Dict, build_dir: str = "../build"):
        """
        Initialize with solver configuration and build directory.
        config: Dictionary containing solver parameters
        build_dir: Path to build directory containing bin/heat_equation
        """
        self.config = config
        self.build_dir = build_dir
        self.executable_path = os.path.join(build_dir, "bin", "heat_equation")

    def validate_environment(self) -> None:
        """Validate that build directory and executable exist."""
        if not os.path.exists(self.build_dir):
            raise FileNotFoundError(
                f"Build directory not found: {self.build_dir}\n"
                "Run 'mkdir build && cd build && cmake .. && make' first"
            )

        if not os.path.exists(self.executable_path):
            raise FileNotFoundError(
                f"Solver executable not found: {self.executable_path}\n"
                "Ensure build is complete with 'make heat_equation'"
            )

    def build_command(self) -> List[str]:
        """
        Construct full mpirun command with all parameters.
        Returns: List of command arguments for subprocess
        """
        # Start with mpirun command
        processes = self.config.get('processes', 1)
        if processes > 1:
            cmd = ['mpirun', '-n', str(processes), self.executable_path]
        else:
            cmd = [self.executable_path]

        # Add solver parameters
        for config_key, cli_arg in self.PARAM_MAPPING.items():
            if config_key in self.config:
                value = self.config[config_key]
                cmd.extend([cli_arg, str(value)])

        return cmd

    def run(self) -> subprocess.Popen:
        """Execute solver and return process handle."""
        self.validate_environment()

        cmd = self.build_command()
        print(f"Running: {' '.join(cmd)}")

        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        return process

    def wait_for_completion(self, process: subprocess.Popen, show_output: bool = False) -> int:
        """
        Wait for solver to finish and optionally show output.
        Returns exit code.
        """
        stdout, stderr = process.communicate()

        if show_output and stdout:
            print(stdout)

        if stderr:
            print(f"Solver stderr: {stderr}")

        return process.returncode

    def get_output_prefix(self) -> str:
        """Get the output file prefix from config."""
        output_prefix = self.config.get('output_prefix', 'output/solution')

        # Ensure output directory exists
        output_dir = os.path.dirname(output_prefix)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

        return output_prefix

    def collect_output_files(self) -> List[str]:
        """
        Find all generated output files matching prefix.
        Returns sorted list of file paths.
        """
        prefix = self.get_output_prefix()

        # Look for files with prefix and .txt or .vtk extensions
        patterns = [
            f"{prefix}*.txt",
            f"{prefix}*.vtk",
            f"{prefix}*.out"
        ]

        files = []
        for pattern in patterns:
            files.extend(glob.glob(pattern))

        return sorted(files)

    def run_and_collect(self, show_output: bool = False) -> List[str]:
        """
        Run solver and collect output files in one call.
        Returns list of generated output files.
        """
        # Get initial file count
        initial_files = set(self.collect_output_files())

        # Run solver
        process = self.run()
        exit_code = self.wait_for_completion(process, show_output)

        if exit_code != 0:
            raise RuntimeError(
                f"Solver exited with code {exit_code}. "
                "Check stderr for details."
            )

        # Collect new files
        final_files = set(self.collect_output_files())
        new_files = list(final_files - initial_files)

        if not new_files:
            print("Warning: No output files found.")
            return []

        print(f"Found {len(new_files)} output file(s)")
        return new_files
```

- [ ] **Step 2: Commit**

```bash
git add visualization/core/solver_runner.py
git commit -m "feat: add MPI solver runner with subprocess execution"
```

---

### Task 7: Create renderers package __init__.py

**Files:**
- Create: `visualization/renderers/__init__.py`

- [ ] **Step 1: Create empty __init__.py file**

Create `visualization/renderers/__init__.py`:
```python
"""
Visualization renderers for heat equation solutions.
"""
```

- [ ] **Step 2: Commit**

```bash
git add visualization/renderers/__init__.py
git commit -m "feat: add renderers package initialization"
```

---

### Task 8: Create heatmap.py rendering module

**Files:**
- Create: `visualization/renderers/heatmap.py`

- [ ] **Step 1: Create heatmap.py with plot_heat_map function**

Create `visualization/renderers/heatmap.py`:
```python
"""
Heat map rendering for single time steps.
"""
import matplotlib.pyplot as plt
import numpy as np
from typing import Optional
from matplotlib.figure import Figure


def plot_heat_map(
    x: np.ndarray,
    y: np.ndarray,
    temperature: np.ndarray,
    title: str = "Heat Distribution",
    colormap: str = "hot",
    show: bool = True,
    save_path: Optional[str] = None,
    colorbar: bool = True,
    dpi: int = 100
) -> Figure:
    """
    Plot 2D heat distribution as pcolormesh.
    Shows colorbar by default.

    Args:
        x: 2D array of x coordinates
        y: 2D array of y coordinates
        temperature: 2D array of temperature values
        title: Plot title
        colormap: Matplotlib colormap name
        show: Display the plot
        save_path: Path to save the figure
        colorbar: Show colorbar
        dpi: Resolution for saved image

    Returns:
        matplotlib Figure for further customization
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # Create heat map
    im = ax.pcolormesh(x, y, temperature, cmap=colormap, shading='auto')

    if colorbar:
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label('Temperature', rotation=270, labelpad=15)

    ax.set_xlabel('x')
    ax.set_ylabel('ℽ')
    ax.set_title(title)
    ax.set_aspect('equal')

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
        print(f"Saved plot to {save_path}")

    if show:
        plt.show()

    return fig
```

- [ ] **Step 2: Commit**

```bash
git add visualization/renderers/heatmap.py
git commit -m "feat: add heat map rendering function"
```

---

### Task 9: Create animation.py module

**Files:**
- Create: `visualization/renderers/animation.py`

- [ ] **Step 1: Create animation.py with create_animation function**

Create `visualization/renderers/animation.py`:
```python
"""
Animation rendering for multiple time steps.
"""
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from typing import List, Dict, Optional
from matplotlib.animation import FuncAnimation


def create_animation(
    files_data: List[Dict],
    output_path: Optional[str] = None,
    show: bool = True,
    fps: int = 10,
    colormap: str = "hot",
    dpi: int = 100
) -> FuncAnimation:
    """
    Create animation from multiple time steps.

    Args:
        files_data: List of parsed solution data from SolutionReader
        output_path: Path to save animation (e.g., 'output.gif')
        show: Display the animation
        fps: Frames per second
        colormap: Matplotlib colormap name
        dpi: Resolution for saved animation

    Returns:
        matplotlib FuncAnimation object
    """
    if not files_data:
        raise ValueError("No data provided for animation")

    # Get grid dimensions from first frame
    nx = files_data[0]['nx']
    ny = files_data[0]['ny']

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Initialize plot with first frame
    x = files_data[0]['x']
    y = files_data[0]['y']
    temperature = files_data[0]['temperature']

    im = ax.pcolormesh(x, y, temperature, cmap=colormap, shading='auto')
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label('Temperature', rotation=270, labelpad=15)

    ax.set_xlabel('x')
    ax.set_ylabel('ℽ')
    ax.set_title(f"Time = {files_data[0]['time']:.4f}, Step = {files_data[0]['step']}")
    ax.set_aspect('equal')

    # Update function for animation
    def update(frame_idx):
        data = files_data[frame_idx]
        im.set_array(data['temperature'].ravel())
        ax.set_title(f"Time = {data['time']:.4f}, Step = {data['step']}")
        return [im]

    # Create animation
    anim = FuncAnimation(
        fig, update, frames=len(files_data),
        interval=1000/fps, blit=True
    )

    if output_path:
        if output_path.endswith('.gif'):
            anim.save(output_path, writer='pillow', fps=fps, dpi=dpi)
        elif output_path.endswith('.mp4'):
            anim.save(output_path, writer='ffmpeg', fps=fps, dpi=dpi)
        else:
            anim.save(output_path + '.gif', writer='pillow', fps=fps, dpi=dpi)
            output_path += '.gif'
        print(f"Saved animation to {output_path}")

    if show:
        plt.show()

    return anim
```

- [ ] **Step 2: Commit**

```bash
git add visualization/renderers/animation.py
git commit -m "feat: add animation creation function"
```

---

### Task 10: Create main entry point run_and_visualize.py

**Files:**
- Create: `visualization/run_and_visualize.py`

- [ ] **Step 1: Create run_and_visualize.py with CLI and workflow orchestration**

Create `visualization/run_and_visualize.py`:
```python
#!/usr/bin/env python3
"""
Run MPI heat equation solver and visualize results automatically.

Usage:
    python run_and_visualize.py --preset medium
    python run_and_visualize.py --preset small --nt 500 --solver CG
"""
import argparse
import yaml
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from core.solver_runner import SolverRunner
from core.solution_reader import SolutionReader
from renderers.heatmap import plot_heat_map
from renderers.animation import create_animation


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run MPI heat equation solver and visualize results",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Preset selection
    parser.add_argument('--preset', choices=['small', 'medium', 'large'],
                        help='Use preset config (small/medium/large)')

    # Solver parameters
    parser.add_argument('-n', '--processes', type=int,
                        help='Number of MPI processes')
    parser.add_argument('-nx', type=int,
                        help='Grid points in x-direction')
    parser.add_argument('-ny', type=int,
                        help='Grid points in y-direction')
    parser.add_argument('-nt', type=int,
                        help='Number of time steps')
    parser.add_argument('-dt', type=float,
                        help='Time step size')
    parser.add_argument('--solver',
                        choices=['Jacobi', 'SOR', 'CG'],
                        help='Solver type')
    parser.add_argument('--scheme',
                        choices=['ImplicitEuler', 'CrankNicolson'],
                        help='Time scheme')

    # Output options
    parser.add_argument('-o', '--output',
                        help='Output file prefix')
    parser.add_argument('--build-dir', default='../build',
                        help='Build directory (default: ../build)')

    # Visualization options
    parser.add_argument('--no-plot', action='store_true',
                        help='Skip plotting')
    parser.add_argument('--save-only', action='store_true',
                        help='Save images, do not show')
    parser.add_argument('--dpi', type=int, default=100,
                        help='Image resolution for saved files')

    # Other
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose output')

    return parser.parse_args()


def load_preset_config(preset_name: str) -> dict:
    """Load preset configuration from YAML file."""
    config_path = os.path.join(
        os.path.dirname(__file__), 'config', 'presets.yaml'
    )

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    if preset_name not in config['presets']:
        raise ValueError(
            f"Preset '{preset_name}' not found. "
            f"Available: {', '.join(config['presets'].keys())}"
        )

    return config['presets'][preset_name].copy()


def merge_config(preset_config: dict, args) -> dict:
    """Merge preset configuration with command-line arguments."""
    result = preset_config.copy() if preset_config else {}

    # Map args to config keys
    arg_mapping = {
        'processes': 'processes',
        'nx': 'nx',
        'ny': 'ny',
        'nt': 'nt',
        'dt': 'dt',
        'solver': 'solver',
        'scheme': 'scheme',
        'output': 'output_prefix',
    }

    for arg_key, config_key in arg_mapping.items():
        value = getattr(args, arg_key)
        if value is not None:
            result[config_key] = value

    return result


def print_config(config: dict):
    """Print current configuration."""
    print("Configuration:")
    for key, value in sorted(config.items()):
        print(f"  {key}: {value}")


def visualize_results(files: list, config: dict, args):
    """Visualize solver results."""
    if not files:
        print("No files to visualize")
        return

    if args.verbose:
        print(f"Visualizing {len(files)} file(s)")

    # Read first file to determine format
    try:
        data = SolutionReader.read_file(files[0])
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    show = not args.save_only

    if len(files) == 1:
        # Single file - create static plot
        if args.verbose:
            print("Creating heat map for single time step")

        title = f"Heat Distribution (Time = {data['time']:.4f}, Step = {data['step']})"
        save = f"{files[0]}.png" if args.save_only else None

        plot_heat_map(
            data['x'], data['y'], data['temperature'],
            title=title,
            show=show,
            save_path=save,
            dpi=args.dpi
        )
    else:
        # Multiple files - create animation
        if args.verbose:
            print(f"Creating animation from {len(files)} time steps")

        # Read all files
        files_data = SolutionReader.read_pattern(
            os.path.join(os.path.dirname(files[0]) or '.', os.path.basename(files[0]).replace('_step', '_step*').replace('.txt', '.txt').replace('.vtk', '.vtk'))
        )

        # Create output path for animation
        output_prefix = config.get('output_prefix', 'output/solution')
        save = f"{output_prefix}_animation.gif" if args.save_only else None

        create_animation(
            files_data,
            output_path=save,
            show=show,
            dpi=args.dpi
        )


def run_workflow(config: dict, args):
    """Execute complete workflow: run solver, monitor, visualize."""
    if args.verbose:
        print_config(config)

    # Create and run solver
    runner = SolverRunner(config, build_dir=args.build_dir)

    try:
        files = runner.run_and_collect(show_output=args.verbose)
    except Exception as e:
        print(f"Error running solver: {e}")
        return

    if not files:
        print("No output files generated")
        return

    # Visualize results
    if not args.no_plot:
        visualize_results(files, config, args)


def main():
    """Main entry point."""
    args = parse_arguments()

    # Load preset configuration if specified
    preset_config = None
    if args.preset:
        preset_config = load_preset_config(args.preset)
        if args.verbose:
            print(f"Loaded preset: {args.preset}")

    # Merge configuration
    config = merge_config(preset_config, args)

    # Run workflow
    run_workflow(config, args)


if __name__ == '__main__':
    main()
```

- [ ] **Step 2: Make script executable**

```bash
chmod +x visualization/run_and_visualize.py
```

- [ ] **Step 3: Commit**

```bash
git add visualization/run_and_visualize.py
git commit -m "feat: add main entry point with CLI and workflow"
```

---

### Task 11: Create test data for manual testing

**Files:**
- Create: `visualization/test_data/test_output.txt`

- [ ] **Step 1: Create sample test output file**

```bash
mkdir -p visualization/test_data
```

Create `visualization/test_data/test_output.txt`:
```txt
# Solution at time = 1.0, step = 100
# Grid: 10 x 10
# Physical domain: [0, 1.0] x [0, 1.0]
# Format: I J x y u(x,y)
0 0 0.0 0.0 0.0
1 0 0.1111111111 0.0 0.05
2 0 0.2222222222 0.0 0.1
3 0 0.3333333333 0.0 0.15
4 0 0.4444444444 0.0 0.2
5 0 0.5555555556 0.0 0.25
6 0 0.6666666667 0.0 0.3
7 0 0.7777777778 0.0 0.35
8 0 0.8888888889 0.0 0.4
9 0 1.0 0.0 0.45
0 1 0.0 0.1111111111 0.05 0.05
1 1 0.1111111111 0.1111111111 0.1
2 1 0.2222222222 0.1111111111 0.15
3 1 0.3333333333 0.1111111111 0.2
4 1 0.4444444444 0.1111111111 0.25
5 1 0.5555555556 0.1111111111 0.3
6 1 0.6666666667 0.1111111111 0.35
7 1 0.7777777778 0.1111111111 0.4
8 1 0.8888888889 0.1111111111 0.45
9 1 1.0 0.1111111111 0.5
0 2 0.0 0.2222222222 0.1 0.1
1 2 0.1111111111 0.2222222222 0.15
2 2 0.2222222222 0.2222222222 0.2
3 2 0.3333333333 0.2222222222 0.25
4 2 0.4444444444 0.2222222222 0.3
5 2 0.5555555556 0.2222222222 0.35
6 2 0.6666666667 0.2222222222 0.4
7 2 0.7777777778 0.2222222222 0.45
8 2 0.8888888889 0.2222222222 0.5
9 2 1.0 0.2222222222 0.55
0 3 0.0 0.3333333333 0.15 0.15
1 3 0.1111111111 0.3333333333 0.2
2 3 0.2222222222 0.3333333333 0.25
3 3 0.3333333333 0.3333333333 0.3
4 3 0.4444444444 0.3333333333 0.35
5 3 0.5555555556 0.3333333333 0.4
6 3x 0.6666666667 0.3333333333 0.45
7 3 0.7777777778 0.3333333333 0.5
8 3 0.8888888889 0.3333333333 0.55
9 3 1.0 0.3333333333 0.6
0 4 0.0 0.4444444444 0.2 0.2
1 4 0.1111111111 0.4444444444 0.25
2 4 0.2222222222 0.4444444444 0.3
3 4 0.3333333333 0.4444444444 0.35
4 4 0.4444444444 0.4444444444 0.4
5 4 0.5555555556 0.4444444444 0.45
6 4 0.6666666667 0.4444444444 0.5
7 4 0.7777777778 0.4444444444 0.55
8 4 0.8888888889 0.4444444444 0.6
9 4 1.0 0.4444444444 0.65
0 5 0.0 0.5555555556 0.25 0.25
1 5 0.1111111111 0.5555555556 0.3
2 5 0.2222222222 0.5555555556 0.35
3 5 0.3333333333 0.5555555556 0.4
4 5 0.4444444444 0.5555555556 0.45
5 5 0.5555555556 0.5555555556 0.5
6 5 0.6666666667 0.5555555556 0.55
7 5 0.7777777778 0.5555555556 0.6
8 5 0.8888888889 0.5555555556 0.65
9 5 1.0 0.5555555556 0.7
0 6 0.0 0.6666666667 0.3 0.3
1 6 0.1111111111 0.6666666667 0.35
2 6 0.2222222222 0.6666666667 0.4
3 6 0.3333333333 0.6666666667 0.45
4 6 0.4444444444 0.6666666667 0.5
5 6 0.5555555556 0.6666666667 0.55
6 6 0.6666666667 0.6666666667 0.6
7 6 0.7777777778 0.6666666667 0.65
8 6 0.8888888889 0.6666666667 0.7
9 6 1.0 0.6666666667 0.75
0 7 0.0 0.7777777778 0.35 0.35
1 7 0.1111111111 0.7777777778 0.4
2 7 0.2222222222 0.7777777778 0.45
3 7 0.3333333333 0.7777777778 0.5
4 7 0.4444444444 0.7777777778 0.55
5 7 0.5555555556 0.7777777778 0.6
6 7 0.6666666667 0.7777777778 0.65
7 7 0.7777777778 0.7777777778 0.7
8 7 0.8888888889 0.7777777778 0.75
9 7 1.0 0.7777777778 0.8
0 8 0.0 0.8888888889 0.4 0.4
1 8 0.1111111111 0.8888888889 0.45
2 8 0.2222222222 0.8888888889 0.5
3 8 0.3333333333 0.8888888889 0.55
4 8 0.4444444444 0.8888888889 0.6
5 8 0.5555555556 0.8888888889 0.65
6 8 0.6666666667 0.8888888889 0.7
7 8 0.7777777778 0.8888888889 0.75
8 8 0.8888888889 0.8888888889 0.8
9 8 1.0 0.8888888889 0.85
0 9 0.0 1.0 0.45 0.45
1 9 0.1111111111 1.0 0.5 0.5
2 9 0.2222222222 1.0 0.55 0.55
3 9 0.3333333333 1.0 0.6 0.6
4 9 0.4444444444 1.0 0.65 0.65
5 9 0.5555555556 1.0 0.7 0.7
6 9 0.6666666667 1.0 0.75 0.75
7 9 0.7777777778 1.0 0.8 0.8
8 9 0.8888888889 1.0 0.85 0.85
9 9 1.0 1.0 0.9 0.9
```

- [ ] **Step 2: Commit**

```bash
git add visualization/test_data/test_output.txt
git commit -m "test: add sample test data for manual testing"
```

---

### Task 12: Create README for visualization

**Files:**
- Create: `visualization/README.md`

- [ ] **Step 1: Create README with usage instructions**

Create `visualization/README.md`:
```markdown
# 2x Heat Equation Solver - Python Visualization

Automated visualization system for the MPI heat equation solver.

## Installation

```bash
cd visualization
pip install -r requirements.txt
```

## Quick Start

### Using Preset Configurations

```bash
# Run with small preset
python run_and_visualize.py --preset small

# Run with medium preset
python run_and_visualize.py --preset medium

# Run with large preset
python run_and_visualize.py --preset large
```

### Custom Parameters

```bash
# Run with custom parameters
python run_and_visualize.py \
    -n 4 \
    -nx 100 \
    -ny 100 \
    -nt 500 \
    -dt 0.005 \
    -solver SOR \
    -scheme CrankNicolson \
    -output results/my_run

# Override preset values
python run_and_visualize.py --preset medium --nt 1000 --solver CG
```

## Command-Line Options

| Option | Description |
|--------|-------------|
| `--preset` | Use preset config (small/medium/large) |
| `-n, --processes` | Number of MPI processes |
| `-nx, -ny` | Grid points in x/y direction |
| `-nt` | Number of time steps |
| `-dt` | Time step size |
| `--solver` | Solver type (Jacobi/SOR/CG) |
| `--scheme` | Time scheme (ImplicitEuler/CrankNicolson) |
| `-o, --output` | Output file prefix |
| `--build-dir` | Build directory (default: ../build) |
| `--no-plot` | Skip plotting |
| `--save-only` | Save images, do not show |
| `--dpi` | Image resolution for saved files |
| `-v, --verbose` | Verbose output |

## Examples

### View Single Time Step

```bash
python run_and_visualize.py --preset small --nt 10
```

### Create Animation

```bash
# Generate multiple outputs (requires solver to output at intervals)
python run_and_visualize.py \
    --preset medium \
    -nt 500 \
    -interval 10  # Note: interval requires solver support
```

### Save Without Display

```bash
python run_and_visualize.py --preset small --save-only --dpi 300
```

## File Formats

The visualization system automatically detects and reads:

- **Text format** (.txt) - Default solver output
- **VTK format** (.vtk) - For ParaView compatibility

## Requirements

- Python 3.7+
- matplotlib >= 3.5.0
- numpy >= 1.21.0
- pyyaml >= 6.0
- tqdm >= 4.62.0
- MPI (for solver execution)

## Testing

To test the visualization with sample data:

```bash
cd visualization
python3 -c "
from core.solution_reader import SolutionReader
from renderers.heatmap import plot_heat_map

data = SolutionReader.read_file('test_data/test_output.txt')
plot_heat_map(data['x'], data['y'], data['temperature'], show=False, save_path='test_plot.png')
print('Test successful!')
"
```
```

- [ ] **Step 2: Commit**

```bash
git add visualization/README.md
git commit -m "docs: add visualization README with usage instructions"
```

---

### Task 13: Update main project README

**Files:**
- Modify: `README.md:456-457` (Add visualization section)

- [ ] **Step 1: Add visualization section to main README**

Add to `README.md` (at the end, before final line):
```markdown
# Visualization

The project includes Python tools for automatic visualization of solver results.

## Quick Visualization

```bash
cd visualization
pip install -r requirements.txt

# Run solver with preset and visualize automatically
python run_and_visualize.py --preset medium

# Custom run with visualization
python run_and_visualize.py -n 4 -nx 100 -ny 100 -nt 500 -solver SOR
```

See [visualization/README.md](visualization/README.md) for detailed usage.
```

- [ ] **Step 2: Commit**

```bash
git add README.md
git commit -m "docs: add visualization section to main README"
```
