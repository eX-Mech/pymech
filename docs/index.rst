.. _documentation-index:

pymech documentation
====================

This is the documentation for pymech_, a Python suite of routines for Nek5000_ and SIMSON_.


+--------------+------------------------------------------------------+
| Authors      |  |author| (:ref:`see here <thanks.md#credits>`)      |
+--------------+------------------------------------------------------+
| Version      |  |version|                                           |
+--------------+------------------------------------------------------+

Pymech can be used for reading, editing and writing Nek5000_ and SIMSON_ meshes
and output files. For a detailed tutorial refer to :ref:`usage
<Usage.ipynb#usage>`.

The data structure is defined by the :py:class:`pymech.exadata.exadata` class, found in :ref:`exadata`.
The functions for manipulating Nek5000_ files are in :ref:`neksuite`, while the
functions for SIMSON_ are, of course, in :ref:`simsonsuite`.

.. _installation:

Installation
------------

Pymech requires Python version 3.7 and above. For most purposes, we recommend
creating a `virtual environment`_ and then running::

   pip install pymech

Optional dependencies can be installed as follows::

   pip install "pymech[full]"

.. _virtual environment: https://packaging.python.org/tutorials/installing-packages/#creating-virtual-environments

-------------------------------------------------------------------------------

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   core.rst
   neksuite.rst
   simsonsuite.rst
   vtksuite.rst
   dataset
   usage.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Reference

   internals
   tests.rst
   deprecations.rst
   changelog.md
   contributing.md
   thanks.md
   release.rst

.. Indices and tables

   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`


.. External links:

.. _pymech: https://github.com/eX-Mech/pymech
.. _SIMSON: https://github.com/KTH-Nek5000/SIMSON
.. _Nek5000: https://nek5000.mcs.anl.gov/
