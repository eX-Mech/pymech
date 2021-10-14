"""The pymech API

.. autosummary::
   :toctree:

   core
   neksuite
   simsonsuite
   vtksuite
   dataset
   log

"""

from .neksuite import *  # noqa
from .simsonsuite import *  # noqa
from .dataset import *  # noqa

from ._version import __version__  # noqa
