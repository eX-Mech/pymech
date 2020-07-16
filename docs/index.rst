.. _documentation-index:

pymech documentation
====================

This is the documentation for pymech_, a Python suite of routines for *Nek5000* and *Simson*.

+--------------+------------------------------+
| Authors      |  | |author|                  |
+--------------+------------------------------+
| Version      |  |release|                   |
+--------------+------------------------------+
| Date         |  |today|                     |
+--------------+------------------------------+
| Installation |   ``pip install pymech``     |
+--------------+------------------------------+

Pymech can be used for reading, editing and writing *Nek5000* and *Simson* output files.

The data structure is defined by the ``exadata`` class, found in :ref:`exadata`.
The functions for manipulating *Nek5000* files are in :ref:`neksuite`, while the
functions for *Simson* are, of course, in :ref:`simsonsuite`.


-------------------------------------------------------------------------------

Contents:
---------

.. toctree::
   :maxdepth: 2

   exadata.rst
   neksuite.rst
   simsonsuite.rst
   vtksuite.rst
   dataset
   tests.rst
   usage.ipynb
   changelog.md



.. Indices and tables

   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`


.. External links:

.. _pymech: https://github.com/jcanton/pymech
