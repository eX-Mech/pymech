.. _dataset:

pymech.dataset
==============

Interface for reading files as xarray_ datasets.

.. _xarray: https://xarray.pydata.org

Installing ``pymech`` also registers as a ``xarray`` backend. This means
in addition to :func:`pymech.dataset.open_dataset`, a file can be directly opened via
:func:`xarray.open_dataset` as follows:

>>> import xarray as xr
>>> xr.open_dataset("case0.f00001")  # let xarray choose the backend / engine
>>> xr.open_dataset("case0.f00001", engine="pymech")  # or explicitly mention the *engine*

.. seealso::

    `The backend API of xarray
    <https://xarray.pydata.org/en/stable/internals/how-to-add-new-backend.html>`__
    and the implementation :class:`PymechXarrayBackend` (:ref:`internals`)

.. todo::

    Opening as a object is not supported by :func:`pymech.neksuite.readnek`


The module also provides :func:`pymech.dataset.open_mfdataset` to open multiple
files and merge into a single dataset. This is a wrapper around
:func:`xarray.open_mfdataset`, with sane defaults such as merge along ``time``
dimension and combine using :func:`xarray.combine_nested`.


Contents of dataset.py
----------------------

.. automodule:: pymech.dataset


.. note::

   See usage_ for more details.

.. _usage: usage.html#pymech-dataset
