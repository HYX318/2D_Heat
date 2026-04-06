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
                        grid_info = parts[1].strip()
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
