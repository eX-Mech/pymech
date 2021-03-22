"""Defines the pymech logger (variable :code:`logger`)."""

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


logger.setLevel(logging.INFO)
