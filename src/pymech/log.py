""":mod:`pymech.log`
=======================

Defines the pymech logger (variable :code:`logger`). The logging level can be
changed to debug level by setting::

    export PYMECH_DEBUG=true

before importing ``pymech``.

For coloured logging::

    pip install rich

"""

import logging
import os
from typing import Union

logger = logging.getLogger("pymech")

# Initialize logging
try:
    # Set a nice colored output
    from rich.logging import RichHandler

    handler: Union[RichHandler, logging.StreamHandler] = RichHandler()
except ImportError:
    # No color available, use default config
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(levelname)s: %(message)s")
    handler.setFormatter(formatter)
    logger.info("Disabling coloured logs; if needed you should `pip install rich`.")

logger.addHandler(handler)

if bool(os.getenv("PYMECH_DEBUG")):
    logger.setLevel(logging.DEBUG)
else:
    logger.setLevel(logging.INFO)
