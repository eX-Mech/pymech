""":mod:`pymech.log`
=======================

Defines the pymech logger (variable :code:`logger`). The logging level can be
changed to debug level by setting::

    export PYMECH_DEBUG=true

before importing ``pymech``.

For coloured logging::

    pip install rich

"""
import os
import logging

logger = logging.getLogger("pymech")

# Initialize logging
try:
    # Set a nice colored output
    from rich.logging import RichHandler

    handler = RichHandler()
    logger.addHandler(handler)
except ImportError:
    # No color available, use default config
    logging.basicConfig(format="%(levelname)s: %(message)s")
    logger.info("Disabling color, you really want to install colorlog.")


if bool(os.getenv("PYMECH_DEBUG")):
    logger.setLevel(logging.DEBUG)
else:
    logger.setLevel(logging.INFO)
