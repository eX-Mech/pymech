"""The pymech API

.. autosummary::
   :toctree:

   core
   neksuite
   simsonsuite
   vtksuite
   dataset
   meshtools
   log
   pyvista_backend

"""

from .neksuite import *  # noqa
from .simsonsuite import *  # noqa

try:
    from .dataset import *  # noqa
except Exception as err:
    import traceback
    from warnings import warn

    traceback.print_exc()
    warn(repr(err), ImportWarning)

from ._version import __version__  # noqa

# Optional visualization backend (PyVista + Matplotlib)
try:
    from . import pyvista_backend  # noqa
except ImportError:
    # PyVista/Matplotlib not installed, visualization features unavailable
    pass
