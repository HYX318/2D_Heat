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
    ax.set_ylabel('y')
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
