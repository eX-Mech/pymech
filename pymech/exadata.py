"""Dummy module for compatibility with pymech version < 1.5.0


.. class:: pymech.exadata.datalims

   Use :class:`pymech.core.DataLims` instead

.. class:: pymech.exadata.elem

   Use :class:`pymech.core.Elem` instead

.. class:: pymech.exadata.exadata

   Use :class:`pymech.core.HexaData` instead

"""
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
