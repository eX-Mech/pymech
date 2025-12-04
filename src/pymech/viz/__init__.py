"""Mesh visualization subpackage with multiple backend support.

This subpackage provides a unified interface for mesh visualization using either
PyVista or Matplotlib backends, with automatic backend selection.

Examples
--------
PyVista backend (interactive 3D):

>>> import pymech as pm
>>> from pymech.viz import plot_mesh
>>>
>>> field = pm.readnek("channel3D_0.f00001")
>>> plot_mesh(field, backend='pyvista')

Matplotlib backend (publication figures):

>>> plot_mesh(field, backend='matplotlib', view='xy')

Auto-select best available backend:

>>> plot_mesh(field, backend='auto')
"""

# Import main API from pyvista_backend module (which is the dispatcher)
from .pyvista_backend import (
    add_boundary_conditions,
    get_available_backends,
    hexa_to_pyvista,
    plot_mesh,
)

# Import Protocol and backend classes for advanced use
from .viz_protocol import DEFAULT_BC_COLORS, MeshBackend

# Optional imports for backend implementations
try:
    from .pyvista_backend_impl import PyVistaBackend
except ImportError:
    PyVistaBackend = None

try:
    from .matplotlib_backend import MatplotlibBackend
except ImportError:
    MatplotlibBackend = None

__all__ = (
    # Main API
    "plot_mesh",
    "get_available_backends",
    # PyVista-specific functions
    "hexa_to_pyvista",
    "add_boundary_conditions",
    # Protocol and colors
    "MeshBackend",
    "DEFAULT_BC_COLORS",
    # Backend classes
    "PyVistaBackend",
    "MatplotlibBackend",
)
