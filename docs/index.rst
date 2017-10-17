====================
pymech documentation
====================

This is the documentation for pymech_, a Python suite of routines for *Nek5000* and *Simson*.

:Authors:
   Jacopo Canton,
   Nicolò Fabbiane

:Version: 1.1 :: 2017/10/17

Pymech can be used for reading, editing and writing *Nek5000* and *Simson* output files.

The data structure is defined by the ``exadata`` class, found in exadata.py_.
The functions for manipulating *Nek5000* files are in neksuite.py_, while the
functions for *Simson* are, of course, in simsonsuite.py_.


-------------------------------------------------------------------------------

Contents:
---------

.. toctree::
   :maxdepth: 2

   exadata.rst
   neksuite.rst
   simsonsuite.rst
   vtksuite.rst
   tests.rst


.. Indices and tables
   
   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`


.. External links:

.. _pymech: https://github.com/jcanton/pymech
.. _exadata.py: exadata.rst
.. _neksuite.py: neksuite.rst
.. _simsonsuite.py: simsonsuite.rst
