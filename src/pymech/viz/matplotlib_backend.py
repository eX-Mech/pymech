"""Matplotlib backend implementation for mesh visualization.

This module provides Matplotlib-based mesh visualization for publication-quality
static figures and basic 3D plots.
"""

from typing import Optional, Tuple, Any
import numpy as np

from .core import HexaData
from .log import logger
from .viz_protocol import (
    MeshBackend,
    DEFAULT_BC_COLORS,
    compute_face_center,
    Resolution,
    View,
    Color,
    Colormap,
)

# Try importing Matplotlib
try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    plt = None
    Axes3D = None

__all__ = ("MatplotlibBackend",)


class MatplotlibBackend:
    """Matplotlib visualization backend implementation."""

    def is_available(self) -> bool:
        """Check if Matplotlib is available."""
        return MATPLOTLIB_AVAILABLE

    def get_backend_name(self) -> str:
        """Get backend name."""
        return "matplotlib"

    def get_capabilities(self) -> dict:
        """Get backend capabilities."""
        return {
            "interactive": False,
            "3d": True,  # Limited 3D support
            "jupyter": True,
            "headless": True,
            "formats": ["png", "jpg", "pdf", "svg", "eps", "ps"],
        }

    def plot_mesh(
        self,
        field: HexaData,
        resolution: Resolution = "linear",
        show_bcs: bool = True,
        bc_field: int = 0,
        show_edges: bool = True,
        style: str = "surface",
        color: Color = None,
        cmap: Colormap = None,
        view: View = None,
        screenshot: Optional[str] = None,
        return_plotter: bool = False,
        figsize: Tuple[float, float] = (10, 8),
        **kwargs,
    ) -> Optional[Any]:
        """Plot mesh using Matplotlib.

        Parameters documented in MeshBackend Protocol.

        Notes
        -----
        Matplotlib backend has limited 3D capabilities compared to PyVista.
        Best used for publication-quality 2D projections.
        """
        if not MATPLOTLIB_AVAILABLE:
            raise ImportError(
                "Matplotlib required for this backend. Install with:\n"
                "    pip install matplotlib"
            )

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')

        # Extract edges from elements
        logger.info("Extracting mesh edges for Matplotlib...")

        for iel, elem in enumerate(field.elem):
            # Get element dimensions
            lx, ly, lz = elem.pos.shape[3], elem.pos.shape[2], elem.pos.shape[1]

            # Define edges based on dimensionality
            if field.ndim == 3:
                edges = _get_hex_edges_3d()
            else:  # 2D
                edges = _get_quad_edges_2d()

            # Get BC for this element if needed
            if show_bcs:
                edge_color = _get_element_bc_color(elem, bc_field)
            elif color:
                edge_color = color
            else:
                edge_color = 'blue'

            # Plot each edge
            for (ix1, iy1, iz1), (ix2, iy2, iz2) in edges:
                p1 = elem.pos[:, iz1, iy1, ix1]
                p2 = elem.pos[:, iz2, iy2, ix2]

                ax.plot(
                    [p1[0], p2[0]],
                    [p1[1], p2[1]],
                    [p1[2], p2[2]],
                    color=edge_color,
                    linewidth=0.5,
                    **kwargs
                )

        # Set labels
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        # Set view angle
        if view == 'xy':
            ax.view_init(elev=90, azim=0)
        elif view == 'xz':
            ax.view_init(elev=0, azim=0)
        elif view == 'yz':
            ax.view_init(elev=0, azim=90)
        elif view:
            # Try to parse as (elev, azim)
            try:
                elev, azim = map(float, view.split(','))
                ax.view_init(elev=elev, azim=azim)
            except (ValueError, AttributeError):
                pass  # Use default view

        # Equal aspect ratio
        _set_axes_equal(ax)

        plt.tight_layout()

        # Save screenshot
        if screenshot:
            plt.savefig(screenshot, dpi=300, bbox_inches='tight')
            logger.info(f"Screenshot saved to {screenshot}")

        if return_plotter:
            return fig
        else:
            plt.show()
            return None


def _get_hex_edges_3d() -> list:
    """Get edge definitions for 3D hexahedron."""
    return [
        # Bottom face
        ((0, 0, 0), (-1, 0, 0)),
        ((-1, 0, 0), (-1, -1, 0)),
        ((-1, -1, 0), (0, -1, 0)),
        ((0, -1, 0), (0, 0, 0)),
        # Top face
        ((0, 0, -1), (-1, 0, -1)),
        ((-1, 0, -1), (-1, -1, -1)),
        ((-1, -1, -1), (0, -1, -1)),
        ((0, -1, -1), (0, 0, -1)),
        # Vertical edges
        ((0, 0, 0), (0, 0, -1)),
        ((-1, 0, 0), (-1, 0, -1)),
        ((-1, -1, 0), (-1, -1, -1)),
        ((0, -1, 0), (0, -1, -1)),
    ]


def _get_quad_edges_2d() -> list:
    """Get edge definitions for 2D quadrilateral."""
    return [
        ((0, 0, 0), (-1, 0, 0)),
        ((-1, 0, 0), (-1, -1, 0)),
        ((-1, -1, 0), (0, -1, 0)),
        ((0, -1, 0), (0, 0, 0)),
    ]


def _get_element_bc_color(elem, bc_field: int) -> tuple:
    """Get representative BC color for an element.

    Since matplotlib edge-based rendering doesn't distinguish faces,
    we use the first non-empty BC color found.
    """
    try:
        for iface in range(elem.bcs.shape[1]):
            bc_type = elem.bcs[bc_field, iface][0]
            if bc_type and bc_type != "E":
                color = DEFAULT_BC_COLORS.get(bc_type)
                if color:
                    return color
    except (IndexError, KeyError):
        pass

    # Default to black
    return (0.0, 0.0, 0.0)


def _set_axes_equal(ax: "Axes3D") -> None:
    """Set 3D plot axes to equal scale.

    Parameters
    ----------
    ax : Axes3D
        Matplotlib 3D axes object
    """
    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])

    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))

    ax.set_xlim3d([origin[0] - radius, origin[0] + radius])
    ax.set_ylim3d([origin[1] - radius, origin[1] + radius])
    ax.set_zlim3d([origin[2] - radius, origin[2] + radius])
