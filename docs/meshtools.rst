.. _meshtools:

pymech.meshtools
================

This module contains some functions to generate and modify *Nek5000* meshes.

.. note::
   Those functions can build the proper connectivity of the elements, storing the elements' neighbours as ``"E  "`` boundary conditions like in the ``.rea`` files. This information seem to not be used at all by Nek5000, and is not contained in the binary ``.re2`` files. In some functions it is possible to discard the connectivity information in order to accelerate the execution. This is particularly relevant in ``merge``.


extrude
-------
This function performs a linear extrusion of a 2D mesh in the z direction. The boundary conditions are kept identical on the existing boundaries and can be specified on the new boundaries.

The mandatory inputs are:

- ``mesh``: a :py:class:`pymech.core.HexaData` structure representing the 2D mesh to extrude;
- ``z``: a list of increasing z coordinates at which the mesh will be extruded.

  Additionally, the following arguments can be specified:

- ``bc1`` and ``bc2``: the boundary conditions to apply on the faces at ``z[0]`` and ``z[-1]``. They are specified as a list of strings, one for each field (velocity, temperature, passive scalars). If nothing is specified, this defaults to periodic boundary conditions for all fields.
- ``internal_bcs``: boolean specifying whether the internal connectivity should be built (``"E  "`` boundary conditions between connected elements). The default is ``True``.

The output is a :py:class:`pymech.core.HexaData` structure representing the extruded 3D mesh.


Contents of meshtools.py
------------------------

.. automodule:: pymech.meshtools

