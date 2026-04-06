"""
MPI solver execution and management.
"""
import subprocess
import os
import glob
from typing import Dict, List, Optional
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

    def sanitize_output_path(self, output_prefix: str) -> str:
        """
        Sanitize output path to prevent directory traversal attacks.
        Rejects paths containing '..' to prevent writing outside intended directory.

        Args:
            output_prefix: The output file prefix path

        Returns:
            Absolute path relative to current working directory

        Raises:
            ValueError: If path contains directory traversal attempts
        """
        # Check for directory traversal attempts
        if '..' in output_prefix:
            raise ValueError(
                f"Path traversal detected in output_prefix: {output_prefix}. "
                "Output paths cannot contain '..' for security reasons."
            )

        # Convert to absolute path to ensure it's within the project structure
        abs_path = os.path.abspath(output_prefix)

        return abs_path

    def validate_parameters(self) -> None:
        """
        Validate all solver parameters in the config.
        Ensures numeric parameters are within reasonable bounds and
        string parameters match expected choices.

        Raises:
            ValueError: If any parameter is invalid
        """
        # Validate numeric parameters
        numeric_params = {
            'nx': (1, 100000),      # Minimum 1, maximum 100000
            'ny': (1, 100000),
            'nt': (1, 1000000),
        }

        for param, (min_val, max_val) in numeric_params.items():
            if param in self.config:
                value = self.config[param]
                try:
                    value = float(value)
                except (TypeError, ValueError):
                    raise ValueError(
                        f"Parameter '{param}' must be numeric, got: {value}"
                    )

                if value < min_val or value > max_val:
                    raise ValueError(
                        f"Parameter '{param}' must be between {min_val} and {max_val}, "
                        f"got: {value}"
                    )

        # Validate dt (time step) - must be positive and reasonable
        if 'dt' in self.config:
            dt = self.config['dt']
            try:
                dt = float(dt)
            except (TypeError, ValueError):
                raise ValueError(f"Parameter 'dt' must be numeric, got: {dt}")

            if dt <= 0:
                raise ValueError(f"Parameter 'dt' must be positive, got: {dt}")

            if dt > 1.0:
                raise ValueError(f"Parameter 'dt' should be less than 1.0, got: {dt}")

        # Validate solver type (string enum) - matching C++ solver options
        valid_solvers = ['Jacobi', 'SOR', 'CG', 'ConjugateGradient', 'GaussSeidel']
        if 'solver' in self.config:
            solver = self.config['solver']
            if solver not in valid_solvers:
                raise ValueError(
                    f"Invalid solver '{solver}'. Must be one of: {', '.join(valid_solvers)}"
                )

        # Validate scheme (string enum) - matching C++ solver options
        valid_schemes = ['ImplicitEuler', 'CrankNicolson', 'ExplicitEuler', 'RungeKutta4']
        if 'scheme' in self.config:
            scheme = self.config['scheme']
            if scheme not in valid_schemes:
                raise ValueError(
                    f"Invalid scheme '{scheme}'. Must be one of: {', '.join(valid_schemes)}"
                )

        # Validate processes count
        if 'processes' in self.config:
            processes = self.config['processes']
            try:
                processes = int(processes)
            except (TypeError, ValueError):
                raise ValueError(f"Parameter 'processes' must be integer, got: {processes}")

            if processes < 1:
                raise ValueError(f"Parameter 'processes' must be at least 1, got: {processes}")

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
        # Validate all parameters before building command
        self.validate_parameters()

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

    def wait_for_completion(
        self,
        process: subprocess.Popen,
        show_output: bool = False,
        timeout: Optional[float] = None
    ) -> int:
        """
        Wait for solver to complete and optionally show output.

        Args:
            process: The subprocess.Popen object
            show_output: Whether to print stdout
            timeout: Maximum time to wait in seconds. None means wait indefinitely.

        Returns:
            Exit code

        Raises:
            subprocess.TimeoutExpired: If process doesn't complete within timeout
        """
        try:
            stdout, stderr = process.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()
            raise subprocess.TimeoutExpired(
                f"Solver process timed out after {timeout} seconds",
                process.pid
            )

        if show_output and stdout:
            print(stdout)

        if stderr:
            print(f"Solver stderr: {stderr}")

        return process.returncode

    def get_output_prefix(self) -> str:
        """Get the output file prefix from config with sanitization."""
        output_prefix = self.config.get('output_prefix', 'output/solution')

        # Sanitize path to prevent directory traversal
        output_prefix = self.sanitize_output_path(output_prefix)

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

    def run_and_collect(
        self,
        show_output: bool = False,
        timeout: Optional[float] = None
    ) -> List[str]:
        """
        Run solver and collect output files in one call.

        Args:
            show_output: Whether to print stdout
            timeout: Maximum time to wait in seconds. None means wait indefinitely.

        Returns:
            List of generated output files

        Raises:
            RuntimeError: If solver exits with non-zero code
            subprocess.TimeoutExpired: If process doesn't complete within timeout
        """
        # Get initial file count
        initial_files = set(self.collect_output_files())

        # Run solver
        process = self.run()
        exit_code = self.wait_for_completion(process, show_output, timeout)

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
