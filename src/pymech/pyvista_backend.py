"""PyVista and Matplotlib-based mesh visualization for pymech

This module provides 3D mesh visualization using PyVista (preferred) or Matplotlib,
optimized for Jupyter notebooks but also supporting standalone scripts.

Examples
--------
PyVista backend (interactive 3D):

>>> import pymech as pm
>>> from pymech.pyvista_backend import plot_mesh
>>>
>>> field = pm.readnek("channel3D_0.f00001")
>>> plot_mesh(field, backend='pyvista')

Matplotlib backend (publication figures):

>>> plot_mesh(field, backend='matplotlib', view='xy')

"""

import warnings
from typing import Optional, Tuple, Dict, Any, Union, Literal
import numpy as np

from .core import HexaData, Elem
from .log import logger

# Try importing PyVista
try:
    import pyvista as pv
    PYVISTA_AVAILABLE = True
except ImportError:
    PYVISTA_AVAILABLE = False
    pv = None

# Try importing Matplotlib
try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    plt = None
    Axes3D = None
    Line3DCollection = None

__all__ = ("plot_mesh", "hexa_to_pyvista", "add_boundary_conditions")

# BC color scheme from meshplot.py (lines 57-69)
DEFAULT_BC_COLORS = {
    "": (0.0, 0.0, 0.0),       # Default/empty - black
    "E": (0.0, 0.0, 0.0),       # Element connectivity - black
    "W": (0.0, 0.0, 0.8),       # Wall - dark blue
    "v": (0.3, 0.3, 1.0),       # Velocity BC - light blue
    "O": (0.8, 0.0, 0.0),       # Outflow - dark red
    "o": (1.0, 0.2, 0.2),       # Outflow variant - red
    "ON": (0.8, 0.4, 0.0),      # Open Neumann - dark orange
    "on": (1.0, 0.6, 0.0),      # Open Neumann variant - orange
    "T": (0.0, 0.8, 0.0),       # Temperature BC - dark green
    "t": (0.3, 1.0, 0.3),       # Temperature variant - green
    "I": (0.95, 0.1, 0.6),      # Insulated - magenta
    "P": (0.5, 0.5, 0.5),       # Periodic - gray
}


def hexa_to_pyvista(
    field: HexaData,
    resolution: Literal["linear", "spectral"] = "linear",
    include_fields: bool = True,
) -> "pv.UnstructuredGrid":
    """Convert HexaData to PyVista UnstructuredGrid.

    Parameters
    ----------
    field : HexaData
        Mesh data structure from pymech
    resolution : {'linear', 'spectral'}, default='linear'
        Mesh resolution strategy:
        - 'linear': Use only corner points (8 vertices per hex, fast)
        - 'spectral': Use all GLL points (subdivides each element, accurate)
    include_fields : bool, default=True
        Whether to include velocity, pressure, temperature as point data

    Returns
    -------
    mesh : pv.UnstructuredGrid
        PyVista mesh with optional field data

    Raises
    ------
    ImportError
        If PyVista is not installed
    ValueError
        If resolution is not 'linear' or 'spectral'

    Notes
    -----
    - Curved edges are approximated with straight lines in 'linear' mode
    - 'spectral' mode subdivides each element into (lx-1)*(ly-1)*(lz-1) sub-cells
    - Field data stored as point data arrays: 'velocity', 'pressure', 'temperature'

    """

    if not PYVISTA_AVAILABLE:
        raise ImportError(
            "PyVista is required for mesh visualization. Install with:\n"
            "    pip install pymech[plot]\n"
            "or:\n"
            "    pip install pyvista"
        )

    if resolution == "linear":
        return _hexa_to_pyvista_linear(field, include_fields)
    elif resolution == "spectral":
        return _hexa_to_pyvista_spectral(field, include_fields)
    else:
        raise ValueError(
            f"resolution must be 'linear' or 'spectral', got '{resolution}'"
        )


def _hexa_to_pyvista_linear(field: HexaData, include_fields: bool) -> "pv.UnstructuredGrid":
    """Convert using only corner vertices (fast, approximate)."""

    nel = field.nel
    ndim = field.ndim

    # Determine cell type and vertex indices
    if ndim == 3:
        nvert = 8
        cell_type = pv.CellType.HEXAHEDRON  # VTK_HEXAHEDRON
        # Vertex ordering for VTK hexahedron
        # (0,0,0), (lx,0,0), (lx,ly,0), (0,ly,0), (0,0,lz), (lx,0,lz), (lx,ly,lz), (0,ly,lz)
        vertex_indices = [
            (0, 0, 0), (-1, 0, 0), (-1, -1, 0), (0, -1, 0),  # bottom face
            (0, 0, -1), (-1, 0, -1), (-1, -1, -1), (0, -1, -1),  # top face
        ]
    else:  # 2D
        nvert = 4
        cell_type = pv.CellType.QUAD
        vertex_indices = [(0, 0, 0), (-1, 0, 0), (-1, -1, 0), (0, -1, 0)]

    # Allocate arrays
    total_points = nel * nvert
    points = np.zeros((total_points, 3), dtype=field.elem[0].pos.dtype)

    # Build connectivity array: [nvert, v0, v1, ..., v_{nvert-1}, nvert, ...]
    cells_list = []
    for i in range(nel):
        cell = [nvert] + list(range(i * nvert, (i + 1) * nvert))
        cells_list.extend(cell)
    cells = np.array(cells_list, dtype=np.int64)

    # Extract corner vertices from each element
    for iel, elem in enumerate(field.elem):
        for ivert, (ix, iy, iz) in enumerate(vertex_indices):
            # elem.pos shape: (3, lz, ly, lx)
            points[iel * nvert + ivert] = elem.pos[:, iz, iy, ix]

    # Create UnstructuredGrid
    cell_types = np.full(nel, cell_type, dtype=np.uint8)
    mesh = pv.UnstructuredGrid(cells, cell_types, points)

    # Add field data if requested
    if include_fields:
        # Velocity field (u, v, w)
        if field.var[1] == 3:
            vel = np.zeros((total_points, 3), dtype=field.elem[0].vel.dtype)
            for iel, elem in enumerate(field.elem):
                for ivert, (ix, iy, iz) in enumerate(vertex_indices):
                    vel[iel * nvert + ivert] = elem.vel[:, iz, iy, ix]
            mesh.point_data["velocity"] = vel
            # Add velocity magnitude for easier visualization
            mesh.point_data["velocity_magnitude"] = np.linalg.norm(vel, axis=1)

        # Pressure field
        if field.var[2] == 1:
            pres = np.zeros(total_points, dtype=field.elem[0].pres.dtype)
            for iel, elem in enumerate(field.elem):
                for ivert, (ix, iy, iz) in enumerate(vertex_indices):
                    pres[iel * nvert + ivert] = elem.pres[0, iz, iy, ix]
            mesh.point_data["pressure"] = pres

        # Temperature field
        if field.var[3] == 1:
            temp = np.zeros(total_points, dtype=field.elem[0].temp.dtype)
            for iel, elem in enumerate(field.elem):
                for ivert, (ix, iy, iz) in enumerate(vertex_indices):
                    temp[iel * nvert + ivert] = elem.temp[0, iz, iy, ix]
            mesh.point_data["temperature"] = temp

    # Add element IDs as cell data
    mesh.cell_data["element_id"] = np.arange(nel)

    return mesh


def _hexa_to_pyvista_spectral(field: HexaData, include_fields: bool) -> "pv.UnstructuredGrid":
    """Convert using all GLL points (slow, accurate for curved elements)."""

    nel = field.nel
    ndim = field.ndim
    lx, ly, lz = field.lr1

    # Points per element and cells per element
    nppel = lx * ly * lz
    if ndim == 3:
        ncpel = (lx - 1) * (ly - 1) * (lz - 1)
        nvert = 8
        cell_type = pv.CellType.HEXAHEDRON
    else:
        ncpel = (lx - 1) * (ly - 1)
        nvert = 4
        cell_type = pv.CellType.QUAD

    total_points = nel * nppel
    total_cells = nel * ncpel

    # Allocate arrays
    points = np.zeros((total_points, 3), dtype=field.elem[0].pos.dtype)
    cells_list = []

    # Extract all GLL points
    for iel, elem in enumerate(field.elem):
        for iz in range(lz):
            for iy in range(ly):
                for ix in range(lx):
                    ipt = iel * nppel + ix + iy * lx + iz * lx * ly
                    points[ipt] = elem.pos[:, iz, iy, ix]

    # Build connectivity for sub-cells
    for iel in range(nel):
        base_pt = iel * nppel
        if ndim == 3:
            for iz in range(lz - 1):
                for iy in range(ly - 1):
                    for ix in range(lx - 1):
                        # Hexahedron vertices (VTK ordering)
                        v0 = base_pt + ix + iy * lx + iz * lx * ly
                        v1 = v0 + 1
                        v2 = v0 + lx + 1
                        v3 = v0 + lx
                        v4 = v0 + lx * ly
                        v5 = v4 + 1
                        v6 = v4 + lx + 1
                        v7 = v4 + lx

                        cells_list.extend([8, v0, v1, v2, v3, v4, v5, v6, v7])
        else:  # 2D
            for iy in range(ly - 1):
                for ix in range(lx - 1):
                    v0 = base_pt + ix + iy * lx
                    v1 = v0 + 1
                    v2 = v0 + lx + 1
                    v3 = v0 + lx

                    cells_list.extend([4, v0, v1, v2, v3])

    cells = np.array(cells_list, dtype=np.int64)
    cell_types = np.full(total_cells, cell_type, dtype=np.uint8)
    mesh = pv.UnstructuredGrid(cells, cell_types, points)

    # Add field data
    if include_fields:
        if field.var[1] == 3:
            vel = np.zeros((total_points, 3), dtype=field.elem[0].vel.dtype)
            for iel, elem in enumerate(field.elem):
                for iz in range(lz):
                    for iy in range(ly):
                        for ix in range(lx):
                            ipt = iel * nppel + ix + iy * lx + iz * lx * ly
                            vel[ipt] = elem.vel[:, iz, iy, ix]
            mesh.point_data["velocity"] = vel
            mesh.point_data["velocity_magnitude"] = np.linalg.norm(vel, axis=1)

        if field.var[2] == 1:
            pres = np.zeros(total_points, dtype=field.elem[0].pres.dtype)
            for iel, elem in enumerate(field.elem):
                for iz in range(lz):
                    for iy in range(ly):
                        for ix in range(lx):
                            ipt = iel * nppel + ix + iy * lx + iz * lx * ly
                            pres[ipt] = elem.pres[0, iz, iy, ix]
            mesh.point_data["pressure"] = pres

        if field.var[3] == 1:
            temp = np.zeros(total_points, dtype=field.elem[0].temp.dtype)
            for iel, elem in enumerate(field.elem):
                for iz in range(lz):
                    for iy in range(ly):
                        for ix in range(lx):
                            ipt = iel * nppel + ix + iy * lx + iz * lx * ly
                            temp[ipt] = elem.temp[0, iz, iy, ix]
            mesh.point_data["temperature"] = temp

    return mesh


def add_boundary_conditions(
    mesh: "pv.UnstructuredGrid",
    field: HexaData,
    bc_field: int = 0,
) -> "pv.PolyData":
    """Add boundary condition information to mesh surface.

    Parameters
    ----------
    mesh : pv.UnstructuredGrid
        Mesh from hexa_to_pyvista
    field : HexaData
        Original data with BC information
    bc_field : int, default=0
        Which BC field to visualize (0=velocity, 1=temperature, etc.)

    Returns
    -------
    surface : pv.PolyData
        Surface mesh with BC data added as cell arrays

    Notes
    -----
    Adds cell data arrays:
    - 'bc_type': String array with BC type for each face
    - 'bc_color': RGB color array for visualization (shape: n_cells x 3)

    """

    if not PYVISTA_AVAILABLE:
        raise ImportError("PyVista required for boundary condition visualization")

    # Extract surface of mesh (external faces only)
    surface = mesh.extract_surface()

    # Initialize BC arrays
    n_faces = surface.n_cells
    bc_types = np.empty(n_faces, dtype='<U3')
    bc_colors = np.zeros((n_faces, 3))

    # Get face centers
    face_centers = surface.cell_centers().points

    # Match surface faces to element faces
    for i in range(n_faces):
        face_center = face_centers[i]
        bc_type = _find_bc_for_face(face_center, field, bc_field)
        bc_types[i] = bc_type
        bc_colors[i] = DEFAULT_BC_COLORS.get(bc_type, (0.5, 0.5, 0.5))

    surface.cell_data["bc_type"] = bc_types
    surface.cell_data["bc_color"] = bc_colors

    return surface


def _find_bc_for_face(center: np.ndarray, field: HexaData, bc_field: int) -> str:
    """Find BC type for a face given its center coordinates."""

    # Tolerance for matching face centers
    tol = 1e-4
    ndim = field.ndim
    nfaces = 2 * ndim  # 4 faces in 2D, 6 faces in 3D

    for iel, elem in enumerate(field.elem):
        # Compute face centers for this element
        for iface in range(nfaces):
            face_center = _compute_face_center(elem, iface, ndim)
            dist = np.linalg.norm(center - face_center)

            if dist < tol:
                try:
                    # BC format: (type, param1, param2, float1, float2, ...)
                    bc = elem.bcs[bc_field, iface][0]
                    return bc if bc else ""
                except (IndexError, KeyError):
                    return ""

    return ""  # No match found


def _compute_face_center(elem: Elem, iface: int, ndim: int) -> np.ndarray:
    """Compute center of a face for an element."""

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


def plot_mesh(
    field: HexaData,
    backend: Literal["pyvista", "matplotlib", "auto"] = "auto",
    resolution: Literal["linear", "spectral"] = "linear",
    show_bcs: bool = True,
    bc_field: int = 0,
    show_edges: bool = True,
    style: str = "surface",
    color: Optional[str] = None,
    cmap: Optional[str] = None,
    jupyter_backend: str = "trame",
    view: Optional[str] = None,
    screenshot: Optional[str] = None,
    return_plotter: bool = False,
    figsize: Tuple[float, float] = (10, 8),
    **kwargs,
) -> Optional[Any]:
    """Plot 3D mesh with boundary conditions using PyVista or Matplotlib.

    Parameters
    ----------
    field : HexaData
        Mesh to visualize
    backend : {'pyvista', 'matplotlib', 'auto'}, default='auto'
        Visualization backend. 'auto' prefers PyVista if available
    resolution : {'linear', 'spectral'}, default='linear'
        Mesh resolution: 'linear' (corners only) or 'spectral' (all GLL points)
    show_bcs : bool, default=True
        Whether to color faces by boundary conditions
    bc_field : int, default=0
        Which BC field to visualize (0=velocity, 1=temperature, ...)
    show_edges : bool, default=True
        Whether to show mesh edges
    style : str, default='surface'
        Visualization style: 'surface', 'wireframe', 'points' (PyVista only)
    color : str, optional
        Uniform color if not showing BCs (e.g., 'white', '#3498db')
    cmap : str, optional
        Colormap for scalar fields (e.g., 'viridis', 'coolwarm')
    jupyter_backend : str, default='trame'
        Backend for Jupyter: 'trame' (interactive), 'static', 'ipyvtklink' (PyVista only)
    view : str, optional
        Camera view: 'xy', 'xz', 'yz', 'iso' (PyVista) or similar for Matplotlib
    screenshot : str, optional
        Save screenshot to this filename
    return_plotter : bool, default=False
        Return plotter/figure object for further customization
    figsize : tuple, default=(10, 8)
        Figure size for Matplotlib backend
    **kwargs
        Additional arguments passed to backend's mesh plotting function

    Returns
    -------
    plotter : pv.Plotter or matplotlib.Figure, optional
        Plotter/figure object if return_plotter=True

    Examples
    --------
    Basic usage with PyVista:

    >>> import pymech as pm
    >>> from pymech.pyvista_backend import plot_mesh
    >>> field = pm.readnek("mesh.nek5000")
    >>> plot_mesh(field)

    Using Matplotlib backend for publication figures:

    >>> plot_mesh(field, backend='matplotlib', view='xy', show_bcs=False)

    In Jupyter notebook with custom camera:

    >>> plotter = plot_mesh(field, return_plotter=True)
    >>> plotter.camera_position = 'xy'
    >>> plotter.show()

    """

    # Determine backend
    if backend == "auto":
        if PYVISTA_AVAILABLE:
            backend = "pyvista"
        elif MATPLOTLIB_AVAILABLE:
            backend = "matplotlib"
        else:
            raise ImportError(
                "Neither PyVista nor Matplotlib is available. Install with:\n"
                "    pip install pymech[plot]"
            )

    if backend == "pyvista":
        return _plot_mesh_pyvista(
            field, resolution, show_bcs, bc_field, show_edges, style,
            color, cmap, jupyter_backend, view, screenshot, return_plotter, **kwargs
        )
    elif backend == "matplotlib":
        return _plot_mesh_matplotlib(
            field, show_bcs, bc_field, show_edges, color, view,
            screenshot, return_plotter, figsize, **kwargs
        )
    else:
        raise ValueError(
            f"backend must be 'pyvista', 'matplotlib', or 'auto', got '{backend}'"
        )


def _plot_mesh_pyvista(
    field, resolution, show_bcs, bc_field, show_edges, style,
    color, cmap, jupyter_backend, view, screenshot, return_plotter, **kwargs
):
    """PyVista backend implementation."""

    if not PYVISTA_AVAILABLE:
        raise ImportError(
            "PyVista required for this backend. Install with:\n"
            "    pip install pymech[plot]"
        )

    # Convert mesh
    logger.info(f"Converting HexaData to PyVista mesh (resolution={resolution})...")
    mesh = hexa_to_pyvista(field, resolution=resolution, include_fields=True)

    # Auto-detect Jupyter environment
    try:
        from IPython import get_ipython
        if get_ipython() is not None and 'IPKernelApp' in get_ipython().config:
            in_notebook = True
        else:
            in_notebook = False
    except (ImportError, AttributeError):
        in_notebook = False

    # Setup plotter
    if in_notebook:
        pv.set_jupyter_backend(jupyter_backend)
        plotter = pv.Plotter(notebook=True)
    else:
        plotter = pv.Plotter()

    # Add mesh with BCs
    if show_bcs:
        logger.info("Extracting boundary conditions...")
        surface = add_boundary_conditions(mesh, field, bc_field)

        # Color by BC
        plotter.add_mesh(
            surface,
            scalars="bc_color",
            rgb=True,
            show_edges=show_edges,
            style=style,
            **kwargs,
        )

        # Add legend for BC types
        _add_bc_legend_pyvista(plotter, surface)
    else:
        # Simple mesh without BC coloring
        mesh_kwargs = {"show_edges": show_edges, "style": style}
        if color:
            mesh_kwargs["color"] = color
        if cmap:
            mesh_kwargs["cmap"] = cmap
        mesh_kwargs.update(kwargs)

        plotter.add_mesh(mesh, **mesh_kwargs)

    # Set camera and axes
    plotter.add_axes()
    if view:
        plotter.camera_position = view
    else:
        plotter.camera_position = 'iso'

    # Show or save
    if screenshot:
        plotter.show(screenshot=screenshot, auto_close=False)
        logger.info(f"Screenshot saved to {screenshot}")

    if return_plotter:
        return plotter
    else:
        plotter.show()
        return None


def _add_bc_legend_pyvista(plotter: "pv.Plotter", surface: "pv.PolyData") -> None:
    """Add legend showing BC types and colors for PyVista."""

    unique_bcs = np.unique(surface.cell_data["bc_type"])

    legend_entries = []
    for bc in unique_bcs:
        if bc and bc in DEFAULT_BC_COLORS:
            color = DEFAULT_BC_COLORS[bc]
            legend_entries.append([bc, color])

    if legend_entries:
        plotter.add_legend(legend_entries, bcolor="white", size=(0.15, 0.15))


def _plot_mesh_matplotlib(
    field, show_bcs, bc_field, show_edges, color, view,
    screenshot, return_plotter, figsize, **kwargs
):
    """Matplotlib backend implementation."""

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
        # Get corner points
        lx, ly, lz = elem.pos.shape[3], elem.pos.shape[2], elem.pos.shape[1]

        if field.ndim == 3:
            # Define 12 edges of a hexahedron
            edges = [
                # Bottom face
                ((0,0,0), (-1,0,0)), ((-1,0,0), (-1,-1,0)), ((-1,-1,0), (0,-1,0)), ((0,-1,0), (0,0,0)),
                # Top face
                ((0,0,-1), (-1,0,-1)), ((-1,0,-1), (-1,-1,-1)), ((-1,-1,-1), (0,-1,-1)), ((0,-1,-1), (0,0,-1)),
                # Vertical edges
                ((0,0,0), (0,0,-1)), ((-1,0,0), (-1,0,-1)), ((-1,-1,0), (-1,-1,-1)), ((0,-1,0), (0,-1,-1)),
            ]
        else:  # 2D
            edges = [
                ((0,0,0), (-1,0,0)), ((-1,0,0), (-1,-1,0)),
                ((-1,-1,0), (0,-1,0)), ((0,-1,0), (0,0,0)),
            ]

        # Plot edges
        for (ix1,iy1,iz1), (ix2,iy2,iz2) in edges:
            p1 = elem.pos[:, iz1, iy1, ix1]
            p2 = elem.pos[:, iz2, iy2, ix2]

            # Determine edge color
            if show_bcs:
                edge_color = 'black'  # Simplified for matplotlib
            elif color:
                edge_color = color
            else:
                edge_color = 'blue'

            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]],
                   color=edge_color, linewidth=0.5, **kwargs)

    # Set labels and view
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    if view == 'xy':
        ax.view_init(elev=90, azim=0)
    elif view == 'xz':
        ax.view_init(elev=0, azim=0)
    elif view == 'yz':
        ax.view_init(elev=0, azim=90)

    # Equal aspect ratio
    _set_axes_equal(ax)

    plt.tight_layout()

    if screenshot:
        plt.savefig(screenshot, dpi=300, bbox_inches='tight')
        logger.info(f"Screenshot saved to {screenshot}")

    if return_plotter:
        return fig
    else:
        plt.show()
        return None


def _set_axes_equal(ax: "Axes3D") -> None:
    """Set 3D plot axes to equal scale."""

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
