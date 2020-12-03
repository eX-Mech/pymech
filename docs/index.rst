.. _documentation-index:

pymech documentation
====================

This is the documentation for pymech_, a Python suite of routines for Nek5000_ and SIMSON_.

+--------------+------------------------------+
| Authors      |  | |author|                  |
+--------------+------------------------------+
| Version      |  |release|                   |
+--------------+------------------------------+
| Date         |  |today|                     |
+--------------+------------------------------+
| Installation |   ``pip install pymech``     |
+--------------+------------------------------+

Pymech can be used for reading, editing and writing Nek5000_ and SIMSON_ output files.

The data structure is defined by the :py:class:`pymech.exadata.exadata` class, found in :ref:`exadata`.
The functions for manipulating Nek5000_ files are in :ref:`neksuite`, while the
functions for SIMSON_ are, of course, in :ref:`simsonsuite`.


-------------------------------------------------------------------------------

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   exadata.rst
   neksuite.rst
   simsonsuite.rst
   vtksuite.rst
   dataset
   usage.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Reference

   tests.rst
   changelog.md
   release.rst

.. Indices and tables

   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`


.. External links:

.. _pymech: https://github.com/jcanton/pymech
.. _SIMSON: https://github.com/KTH-Nek5000/SIMSON
.. _Nek5000: https://nek5000.mcs.anl.gov/
