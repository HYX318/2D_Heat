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
