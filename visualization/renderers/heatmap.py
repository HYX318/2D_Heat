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
        show: Display plot
        save_path: Path to save figure
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
    ax.set_ylabel('y')
    ax.set_title(title)
    ax.set_aspect('equal')

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
        print(f"Saved plot to {save_path}")

    if show:
        plt.show()

    return fig
