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
    os.makedirs(directory, exist_ok=True)

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
    os.makedirs(directory, exist_ok=True)

    start_time = time.time()

    while time.time() - start_time < timeout:
        files = glob.glob(os.path.join(directory, pattern))
        if len(files) >= expected_count:
            return sorted(files)
        time.sleep(0.5)

    # Timeout reached
    current_files = glob.glob(os.path.join(directory, pattern))
    return sorted(current_files)
