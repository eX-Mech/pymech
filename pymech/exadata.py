"""Dummy module for compatibility with pymech <= 1.4.1"""
from warnings import warn

from pymech.core import DataLims as datalims  # noqa
from pymech.core import Elem as elem  # noqa
from pymech.core import HexaData as exadata  # noqa


warn(
    (
        "Module pymech.exadata is now pymech.core. This module is kept for backwards "
        "compatibility and would disappear in version 2.0.0"
    ),
    DeprecationWarning,
)
