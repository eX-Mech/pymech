"""The pymech API

.. autosummary::
   :toctree:

   exadata
   neksuite
   simsonsuite
   vtksuite
   dataset
   log

"""

from . import neksuite
from . import simsonsuite

__all__ = ["neksuite", "simsonsuite"]

from ._version import __version__
