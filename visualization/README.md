# 2D Heat Equation Solver - Python Visualization

Automated visualization system for MPI heat equation solver.

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
| `--` | Image resolution for saved files |
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
