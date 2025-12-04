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
