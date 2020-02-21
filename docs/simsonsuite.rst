.. _simsonsuite:

simsonsuite.py
==============

This module contains the functions used to interact with *Simson* files.
The :ref:`simsonsuite-contents` are reported at the bottom of this page.


readdns
-------
This function reads the binary (``.u``) files that *Simson* uses to store
output flow fields.

The only input needed by this function is:

- ``fname``: a string containing the name of the file;

readdns() is clever enough to figure out the rest.
The output is a single ``exadata`` (:ref:`exadata`) structure, with just one
element, containing all the information that was stored in the file.

readplane
---------
This function reads the binary (``.stat``) files that *Simson* uses to store
output statistics file, to be read with ``pxyst``.

The only input needed by this function is:

- ``fname``: a string containing the name of the file;

readplane() is clever enough to figure out the rest.
The output consists of:

========   =================   ===========================================
``x``      array of floats     coordinates of the grid points
``d``      array of floats     u,v,w velocity at the grid points
``nn``     array of integers   number of grid points along each coordinate
``ndim``   integer             number of spatial dimensions of the data
========   =================   ===========================================


-------------------------------------------------------------------------------

.. _simsonsuite-contents:

Contents of simsonsuite.py
--------------------------

**readdns()**

.. literalinclude:: ../pymech/simsonsuite.py
   :language: python
   :lines: 15-241

**readplane()**

.. literalinclude:: ../pymech/simsonsuite.py
   :language: python
   :lines: 245-363
