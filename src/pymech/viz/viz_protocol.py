"""Protocol definition for mesh visualization backends.

This module defines the interface that all visualization backends must implement,
using typing.Protocol for structural subtyping.
"""

from typing import Any, Literal, Optional, Protocol, Tuple, runtime_checkable

import numpy as np
from typing_extensions import TypeAlias

from ..core import HexaData

# Type aliases
Resolution: TypeAlias = Literal["linear", "spectral"]
View: TypeAlias = Optional[str]
Color: TypeAlias = Optional[str]
Colormap: TypeAlias = Optional[str]

# BC color scheme - shared across all backends
DEFAULT_BC_COLORS = {
    "": (0.0, 0.0, 0.0),  # Default/empty - black
    "E": (0.0, 0.0, 0.0),  # Element connectivity - black
    "W": (0.0, 0.0, 0.8),  # Wall - dark blue
    "v": (0.3, 0.3, 1.0),  # Velocity BC - light blue
    "O": (0.8, 0.0, 0.0),  # Outflow - dark red
    "o": (1.0, 0.2, 0.2),  # Outflow variant - red
    "ON": (0.8, 0.4, 0.0),  # Open Neumann - dark orange
    "on": (1.0, 0.6, 0.0),  # Open Neumann variant - orange
    "T": (0.0, 0.8, 0.0),  # Temperature BC - dark green
    "t": (0.3, 1.0, 0.3),  # Temperature variant - green
    "I": (0.95, 0.1, 0.6),  # Insulated - magenta
    "P": (0.5, 0.5, 0.5),  # Periodic - gray
}


@runtime_checkable
class MeshBackend(Protocol):
    """Protocol for mesh visualization backends.

    All visualization backends must implement this interface to ensure
    consistent API across different rendering engines.
    """

    def is_available(self) -> bool:
        """Check if this backend is available (dependencies installed).

        Returns
        -------
        bool
            True if backend can be used, False otherwise
        """
        ...

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
        """Plot mesh with boundary conditions.

        Parameters
        ----------
        field : HexaData
            Mesh to visualize
        resolution : {'linear', 'spectral'}, default='linear'
            Mesh resolution
        show_bcs : bool, default=True
            Whether to color by boundary conditions
        bc_field : int, default=0
            Which BC field to visualize
        show_edges : bool, default=True
            Whether to show mesh edges
        style : str, default='surface'
            Visualization style
        color : str, optional
            Uniform color
        cmap : str, optional
            Colormap for scalar fields
        view : str, optional
            Camera view angle
        screenshot : str, optional
            Save screenshot to file
        return_plotter : bool, default=False
            Return plotter/figure object
        figsize : tuple, default=(10, 8)
            Figure size
        **kwargs
            Backend-specific options

        Returns
        -------
        plotter : optional
            Plotter/figure object if return_plotter=True
        """
        ...

    def get_backend_name(self) -> str:
        """Get the name of this backend.

        Returns
        -------
        str
            Backend name (e.g., 'pyvista', 'matplotlib')
        """
        ...

    def get_capabilities(self) -> dict:
        """Get backend capabilities.

        Returns
        -------
        dict
            Dictionary describing backend features:
            - 'interactive': bool - supports interactive manipulation
            - '3d': bool - supports true 3D rendering
            - 'jupyter': bool - works in Jupyter notebooks
            - 'headless': bool - supports headless rendering
            - 'formats': list - supported export formats
        """
        ...


def get_bc_color(bc_type: str) -> Tuple[float, float, float]:
    """Get RGB color for a boundary condition type.

    Parameters
    ----------
    bc_type : str
        Boundary condition type identifier

    Returns
    -------
    tuple
        RGB color tuple (values 0-1)
    """
    return DEFAULT_BC_COLORS.get(bc_type, (0.5, 0.5, 0.5))


def compute_face_center(elem, iface: int, ndim: int) -> np.ndarray:
    """Compute center of a face for an element.

    This is a shared utility function used by multiple backends.

    Parameters
    ----------
    elem : Elem
        Element object
    iface : int
        Face index
    ndim : int
        Number of dimensions (2 or 3)

    Returns
    -------
    np.ndarray
        3D coordinates of face center
    """
    lx, ly, lz = elem.pos.shape[3], elem.pos.shape[2], elem.pos.shape[1]

    if ndim == 3:
        # 3D: 6 faces (x-, x+, y-, y+, z-, z+)
        if iface == 0:  # x- face
            face_pts = elem.pos[:, :, :, 0]
        elif iface == 1:  # x+ face
            face_pts = elem.pos[:, :, :, -1]
        elif iface == 2:  # y- face
            face_pts = elem.pos[:, :, 0, :]
        elif iface == 3:  # y+ face
            face_pts = elem.pos[:, :, -1, :]
        elif iface == 4:  # z- face
            face_pts = elem.pos[:, 0, :, :]
        else:  # iface == 5, z+ face
            face_pts = elem.pos[:, -1, :, :]
    else:  # 2D
        # 2D: 4 faces (x-, x+, y-, y+)
        if iface == 0:  # x- face
            face_pts = elem.pos[:, 0, :, 0]
        elif iface == 1:  # x+ face
            face_pts = elem.pos[:, 0, :, -1]
        elif iface == 2:  # y- face
            face_pts = elem.pos[:, 0, 0, :]
        else:  # iface == 3, y+ face
            face_pts = elem.pos[:, 0, -1, :]

    # Return mean of all points on the face
    return face_pts.mean(axis=tuple(range(1, face_pts.ndim)))
