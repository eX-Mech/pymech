"""PyVista backend implementation for mesh visualization.

This module provides PyVista-specific mesh visualization, optimized for
interactive 3D rendering in Jupyter notebooks.
"""

from typing import Any, Optional, Tuple

import numpy as np

from ..core import HexaData
from ..log import logger
from .viz_protocol import (
    DEFAULT_BC_COLORS,
    Color,
    Colormap,
    Resolution,
    View,
    compute_face_center,
)

# Try importing PyVista
try:
    import pyvista as pv

    PYVISTA_AVAILABLE = True
except ImportError:
    PYVISTA_AVAILABLE = False
    pv = None

__all__ = ("PyVistaBackend", "hexa_to_pyvista", "add_boundary_conditions")


class PyVistaBackend:
    """PyVista visualization backend implementation."""

    def is_available(self) -> bool:
        """Check if PyVista is available."""
        return PYVISTA_AVAILABLE

    def get_backend_name(self) -> str:
        """Get backend name."""
        return "pyvista"

    def get_capabilities(self) -> dict:
        """Get backend capabilities."""
        return {
            "interactive": True,
            "3d": True,
            "jupyter": True,
            "headless": True,
            "formats": ["png", "jpg", "bmp", "tif", "svg", "eps", "ps", "pdf", "tex"],
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
        """Plot mesh using PyVista.

        Parameters documented in MeshBackend Protocol.
        """
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

            if get_ipython() is not None and "IPKernelApp" in get_ipython().config:
                in_notebook = True
            else:
                in_notebook = False
        except (ImportError, AttributeError):
            in_notebook = False

        # Extract jupyter_backend from kwargs
        jupyter_backend = kwargs.pop("jupyter_backend", "trame")

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
            _add_bc_legend(plotter, surface)
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
            plotter.camera_position = "iso"

        # Show or save
        if screenshot:
            plotter.show(screenshot=screenshot, auto_close=False)
            logger.info(f"Screenshot saved to {screenshot}")

        if return_plotter:
            return plotter
        else:
            plotter.show()
            return None


def hexa_to_pyvista(
    field: HexaData,
    resolution: Resolution = "linear",
    include_fields: bool = True,
) -> "pv.UnstructuredGrid":
    """Convert HexaData to PyVista UnstructuredGrid.

    Parameters
    ----------
    field : HexaData
        Mesh data structure from pymech
    resolution : {'linear', 'spectral'}, default='linear'
        Mesh resolution strategy
    include_fields : bool, default=True
        Whether to include velocity, pressure, temperature as point data

    Returns
    -------
    mesh : pv.UnstructuredGrid
        PyVista mesh with optional field data
    """
    if not PYVISTA_AVAILABLE:
        raise ImportError("PyVista required. Install with: pip install pyvista")

    if resolution == "linear":
        return _hexa_to_pyvista_linear(field, include_fields)
    elif resolution == "spectral":
        return _hexa_to_pyvista_spectral(field, include_fields)
    else:
        raise ValueError(
            f"resolution must be 'linear' or 'spectral', got '{resolution}'"
        )


def _hexa_to_pyvista_linear(
    field: HexaData, include_fields: bool
) -> "pv.UnstructuredGrid":
    """Convert using only corner vertices (fast, approximate)."""
    nel = field.nel
    ndim = field.ndim

    # Determine cell type and vertex indices
    if ndim == 3:
        nvert = 8
        cell_type = pv.CellType.HEXAHEDRON
        vertex_indices = [
            (0, 0, 0),
            (-1, 0, 0),
            (-1, -1, 0),
            (0, -1, 0),  # bottom face
            (0, 0, -1),
            (-1, 0, -1),
            (-1, -1, -1),
            (0, -1, -1),  # top face
        ]
    else:  # 2D
        nvert = 4
        cell_type = pv.CellType.QUAD
        vertex_indices = [(0, 0, 0), (-1, 0, 0), (-1, -1, 0), (0, -1, 0)]

    # Allocate arrays
    total_points = nel * nvert
    points = np.zeros((total_points, 3), dtype=field.elem[0].pos.dtype)

    # Build connectivity
    cells_list = []
    for i in range(nel):
        cell = [nvert] + list(range(i * nvert, (i + 1) * nvert))
        cells_list.extend(cell)
    cells = np.array(cells_list, dtype=np.int64)

    # Extract corner vertices
    for iel, elem in enumerate(field.elem):
        for ivert, (ix, iy, iz) in enumerate(vertex_indices):
            points[iel * nvert + ivert] = elem.pos[:, iz, iy, ix]

    # Create UnstructuredGrid
    cell_types = np.full(nel, cell_type, dtype=np.uint8)
    mesh = pv.UnstructuredGrid(cells, cell_types, points)

    # Add field data
    if include_fields:
        _add_field_data(mesh, field, vertex_indices, nvert, total_points)

    # Add element IDs
    mesh.cell_data["element_id"] = np.arange(nel)

    return mesh


def _hexa_to_pyvista_spectral(
    field: HexaData, include_fields: bool
) -> "pv.UnstructuredGrid":
    """Convert using all GLL points (slow, accurate)."""
    nel = field.nel
    ndim = field.ndim
    lx, ly, lz = field.lr1

    nppel = lx * ly * lz
    if ndim == 3:
        ncpel = (lx - 1) * (ly - 1) * (lz - 1)
        cell_type = pv.CellType.HEXAHEDRON
    else:
        ncpel = (lx - 1) * (ly - 1)
        cell_type = pv.CellType.QUAD

    total_points = nel * nppel
    total_cells = nel * ncpel

    points = np.zeros((total_points, 3), dtype=field.elem[0].pos.dtype)
    cells_list = []

    # Extract all GLL points
    for iel, elem in enumerate(field.elem):
        for iz in range(lz):
            for iy in range(ly):
                for ix in range(lx):
                    ipt = iel * nppel + ix + iy * lx + iz * lx * ly
                    points[ipt] = elem.pos[:, iz, iy, ix]

    # Build connectivity
    for iel in range(nel):
        base_pt = iel * nppel
        if ndim == 3:
            for iz in range(lz - 1):
                for iy in range(ly - 1):
                    for ix in range(lx - 1):
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
        _add_spectral_field_data(mesh, field, lx, ly, lz, nppel, total_points)

    return mesh


def _add_field_data(mesh, field, vertex_indices, nvert, total_points):
    """Add field data to linear mesh."""
    # Velocity
    if field.var[1] == 3:
        vel = np.zeros((total_points, 3), dtype=field.elem[0].vel.dtype)
        for iel, elem in enumerate(field.elem):
            for ivert, (ix, iy, iz) in enumerate(vertex_indices):
                vel[iel * nvert + ivert] = elem.vel[:, iz, iy, ix]
        mesh.point_data["velocity"] = vel
        mesh.point_data["velocity_magnitude"] = np.linalg.norm(vel, axis=1)

    # Pressure
    if field.var[2] == 1:
        pres = np.zeros(total_points, dtype=field.elem[0].pres.dtype)
        for iel, elem in enumerate(field.elem):
            for ivert, (ix, iy, iz) in enumerate(vertex_indices):
                pres[iel * nvert + ivert] = elem.pres[0, iz, iy, ix]
        mesh.point_data["pressure"] = pres

    # Temperature
    if field.var[3] == 1:
        temp = np.zeros(total_points, dtype=field.elem[0].temp.dtype)
        for iel, elem in enumerate(field.elem):
            for ivert, (ix, iy, iz) in enumerate(vertex_indices):
                temp[iel * nvert + ivert] = elem.temp[0, iz, iy, ix]
        mesh.point_data["temperature"] = temp


def _add_spectral_field_data(mesh, field, lx, ly, lz, nppel, total_points):
    """Add field data to spectral mesh."""
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
        Which BC field to visualize

    Returns
    -------
    surface : pv.PolyData
        Surface mesh with BC data
    """
    if not PYVISTA_AVAILABLE:
        raise ImportError("PyVista required")

    # Extract surface
    surface = mesh.extract_surface()

    # Initialize BC arrays
    n_faces = surface.n_cells
    bc_types = np.empty(n_faces, dtype="<U3")
    bc_colors = np.zeros((n_faces, 3))

    # Get face centers
    face_centers = surface.cell_centers().points

    # Match faces to BCs
    for i in range(n_faces):
        face_center = face_centers[i]
        bc_type = _find_bc_for_face(face_center, field, bc_field)
        bc_types[i] = bc_type
        bc_colors[i] = DEFAULT_BC_COLORS.get(bc_type, (0.5, 0.5, 0.5))

    surface.cell_data["bc_type"] = bc_types
    surface.cell_data["bc_color"] = bc_colors

    return surface


def _find_bc_for_face(center: np.ndarray, field: HexaData, bc_field: int) -> str:
    """Find BC type for a face."""
    tol = 1e-4
    ndim = field.ndim
    nfaces = 2 * ndim

    for iel, elem in enumerate(field.elem):
        for iface in range(nfaces):
            face_center = compute_face_center(elem, iface, ndim)
            dist = np.linalg.norm(center - face_center)

            if dist < tol:
                try:
                    bc = elem.bcs[bc_field, iface][0]
                    return bc if bc else ""
                except (IndexError, KeyError):
                    return ""

    return ""


def _add_bc_legend(plotter: "pv.Plotter", surface: "pv.PolyData") -> None:
    """Add BC legend to plotter."""
    unique_bcs = np.unique(surface.cell_data["bc_type"])

    legend_entries = []
    for bc in unique_bcs:
        if bc and bc in DEFAULT_BC_COLORS:
            color = DEFAULT_BC_COLORS[bc]
            legend_entries.append([bc, color])

    if legend_entries:
        plotter.add_legend(legend_entries, bcolor="white", size=(0.15, 0.15))
