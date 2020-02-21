.. _neksuite:

neksuite.py
===========

This module contains the functions used to interact with *Nek5000* files.
The :ref:`neksuite-contents` are reported at the bottom of this page.


readnek
-------
This function reads the binary ``.f%05d`` files that *Nek5000* uses to store
output flow fields.

The only input needed by this function is:

- ``fname``: a string containing the name of the file;

readnek() is clever enough to figure out the rest.
The output is a single ``exadata`` (:ref:`exadata`) structure containing all
the information that was stored in the file.


writenek
--------
This function writes a binary ``.f%05d`` file in the same format that *Nek5000*
uses to store output flow fields. It is therefore possible to open this file with
*VisIt* or *Paraview*.

The inputs needed by this function are:

- ``fname``: a string containing the name of the file;
- ``data``: a single ``exadata`` (:ref:`exadata`) structure containing the
  data to be written to file.

writenek() produces no output.


readrea
-------
This function reads an ASCII ``.rea`` file that *Nek5000* uses to store
simulation parameters and mesh (when a ``.re2`` file is absent).
**ONLY the mesh is actually read**, all parameters are disregarded.

The only input needed by this function is:

- ``fname``: a string containing the name of the file;

readrea() is clever enough to figure out the rest.
The output is a single ``exadata`` (:ref:`exadata`) structure containing the
mesh that was stored in the file.


writerea
--------
This function writes an ASCII ``.rea`` file in the same format that *Nek5000*
uses to store simulation parameters and mesh.
**ONLY the mesh is actually written**, all simulation parameters are initialised
to sensible defaults.

The inputs needed by this function are:

- ``fname``: a string containing the name of the file;
- ``data``: a single ``exadata`` (:ref:`exadata`) structure containing the
  mesh to be written to file.

writerea() produces no output.


-------------------------------------------------------------------------------

.. _neksuite-contents:

Contents of neksuite.py
-----------------------

**readnek()**

.. literalinclude:: ../pymech/neksuite.py
   :language: python
   :lines: 16-207

**writenek()**

.. literalinclude:: ../pymech/neksuite.py
   :language: python
   :lines: 211-348

**readrea()**

.. literalinclude:: ../pymech/neksuite.py
   :language: python
   :lines: 352-512

**writerea()**

.. literalinclude:: ../pymech/neksuite.py
   :language: python
   :lines: 516-770
