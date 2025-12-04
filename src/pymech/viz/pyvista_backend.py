"""Unified mesh visualization API with multiple backend support.

This module provides a unified interface for mesh visualization using either
PyVista or Matplotlib backends, with automatic backend selection.

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

Auto-select best available backend:

>>> plot_mesh(field, backend='auto')

"""

from typing import Optional, Tuple, Any, Literal
import warnings

from .core import HexaData
from .log import logger
from .viz_protocol import MeshBackend, DEFAULT_BC_COLORS

# Import backend implementations
try:
    from .pyvista_backend_impl import PyVistaBackend, hexa_to_pyvista, add_boundary_conditions
    PYVISTA_BACKEND_AVAILABLE = True
except ImportError:
    PYVISTA_BACKEND_AVAILABLE = False
    hexa_to_pyvista = None
    add_boundary_conditions = None

try:
    from .matplotlib_backend import MatplotlibBackend
    MATPLOTLIB_BACKEND_AVAILABLE = True
except ImportError:
    MATPLOTLIB_BACKEND_AVAILABLE = False

__all__ = ("plot_mesh", "hexa_to_pyvista", "add_boundary_conditions", "get_available_backends")


def get_available_backends() -> dict:
    """Get dictionary of available visualization backends.

    Returns
    -------
    dict
        Mapping of backend names to backend instances (only available backends)

    Examples
    --------
    >>> backends = get_available_backends()
    >>> print(f"Available backends: {list(backends.keys())}")
    Available backends: ['pyvista', 'matplotlib']
    """
    backends = {}

    if PYVISTA_BACKEND_AVAILABLE:
        pv_backend = PyVistaBackend()
        if pv_backend.is_available():
            backends['pyvista'] = pv_backend

    if MATPLOTLIB_BACKEND_AVAILABLE:
        mpl_backend = MatplotlibBackend()
        if mpl_backend.is_available():
            backends['matplotlib'] = mpl_backend

    return backends


def _get_backend(backend_name: str) -> MeshBackend:
    """Get a specific backend instance.

    Parameters
    ----------
    backend_name : str
        Name of backend ('pyvista', 'matplotlib', or 'auto')

    Returns
    -------
    MeshBackend
        Backend instance

    Raises
    ------
    ValueError
        If backend name is invalid
    ImportError
        If requested backend is not available
    """
    available_backends = get_available_backends()

    if backend_name == 'auto':
        # Prefer PyVista if available
        if 'pyvista' in available_backends:
            return available_backends['pyvista']
        elif 'matplotlib' in available_backends:
            return available_backends['matplotlib']
        else:
            raise ImportError(
                "No visualization backend available. Install with:\n"
                "    pip install pymech[plot]"
            )
    elif backend_name in available_backends:
        return available_backends[backend_name]
    elif backend_name in ('pyvista', 'matplotlib'):
        # Backend name is valid but not available
        raise ImportError(
            f"{backend_name} backend not available. Install with:\n"
            f"    pip install pymech[plot]"
        )
    else:
        available = list(available_backends.keys())
        raise ValueError(
            f"Invalid backend '{backend_name}'. "
            f"Must be 'pyvista', 'matplotlib', or 'auto'. "
            f"Available: {available}"
        )


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

    This is the main entry point for mesh visualization. It automatically
    selects the best available backend or uses the specified one.

    Parameters
    ----------
    field : HexaData
        Mesh to visualize
    backend : {'pyvista', 'matplotlib', 'auto'}, default='auto'
        Visualization backend. 'auto' prefers PyVista if available
    resolution : {'linear', 'spectral'}, default='linear'
        Mesh resolution: 'linear' (corners only) or 'spectral' (all GLL points)
        Note: Only affects PyVista backend
    show_bcs : bool, default=True
        Whether to color faces/edges by boundary conditions
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
        Backend for Jupyter: 'trame' (interactive), 'static', 'ipyvtklink'
        (PyVista only)
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

    Raises
    ------
    ImportError
        If no visualization backend is available
    ValueError
        If invalid backend name is specified

    Examples
    --------
    Basic usage with auto backend selection:

    >>> import pymech as pm
    >>> from pymech.pyvista_backend import plot_mesh
    >>> field = pm.readnek("mesh.nek5000")
    >>> plot_mesh(field)

    Using specific backend:

    >>> plot_mesh(field, backend='matplotlib', view='xy')

    Customizing PyVista visualization:

    >>> plotter = plot_mesh(field, backend='pyvista', return_plotter=True)
    >>> plotter.camera_position = [(10, 10, 10), (0, 0, 0), (0, 1, 0)]
    >>> plotter.show()

    Saving high-resolution figure:

    >>> plot_mesh(field, backend='matplotlib', screenshot='mesh.pdf',
    ...           figsize=(12, 10))

    See Also
    --------
    get_available_backends : Check which backends are available
    hexa_to_pyvista : Convert HexaData to PyVista format (PyVista backend only)
    add_boundary_conditions : Add BC data to mesh (PyVista backend only)

    """
    # Get backend instance
    backend_instance = _get_backend(backend)

    logger.info(f"Using {backend_instance.get_backend_name()} backend for visualization")

    # Call backend's plot_mesh method
    return backend_instance.plot_mesh(
        field=field,
        resolution=resolution,
        show_bcs=show_bcs,
        bc_field=bc_field,
        show_edges=show_edges,
        style=style,
        color=color,
        cmap=cmap,
        view=view,
        screenshot=screenshot,
        return_plotter=return_plotter,
        figsize=figsize,
        jupyter_backend=jupyter_backend,  # Pass through kwargs
        **kwargs,
    )


# Maintain backward compatibility: export PyVista-specific functions if available
if not PYVISTA_BACKEND_AVAILABLE:
    def hexa_to_pyvista(*args, **kwargs):
        """PyVista not available."""
        raise ImportError(
            "PyVista backend not available. Install with:\n"
            "    pip install pymech[plot]"
        )

    def add_boundary_conditions(*args, **kwargs):
        """PyVista not available."""
        raise ImportError(
            "PyVista backend not available. Install with:\n"
            "    pip install pymech[plot]"
        )


# Module-level convenience: show available backends on import
def _show_backend_info():
    """Display information about available backends (suppressed by default)."""
    backends = get_available_backends()
    if backends:
        logger.debug(f"Available visualization backends: {list(backends.keys())}")
    else:
        warnings.warn(
            "No visualization backends available. Install with: pip install pymech[plot]",
            ImportWarning,
            stacklevel=2
        )


# Don't show info by default to avoid clutter
# _show_backend_info()
