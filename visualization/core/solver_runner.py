"""
MPI solver execution and management.
"""
import subprocess
import os
import glob
from typing import Dict, List
from .output_monitor import monitor_directory


class SolverRunner:
    """Manages execution of MPI heat equation solver."""

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
        Wait for solver to complete and optionally show output.
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
