.. _core:

pymech.core
===========

This module contains the classes for the data structures used by
:ref:`neksuite` and :ref:`simsonsuite`.
The :ref:`core-contents` are reported at the bottom of this page.

HexaData
--------

The main data class is :py:class:`pymech.core.HexaData` which is a structure
for general and element by element information on data stored on hexahedral
element(s).  The general information is constituted by:

==========   ============================   ===============================================
``ndim``     integer                        number of spatial dimensions of the data
``nel``      integer                        number of elements included in :py:class:`pymech.core.HexaData`
``ncurv``    integer                        number of curved sides (used only by *Nek5000*)
``var``      string                         variables included in :py:class:`pymech.core.HexaData` (e.g. 'XUPT')
``lr1``      ``ndim``-array of integers     number of datapoints per element
``time``     float                          simulation time of the file
``istep``    integer                        simulation time step of the file
``wdsz``     integer                        word size, i.e. double or single precision file
``endian``   string                         endianness_ of the file (little/big)
``lims``     :ref:`datalims`                extrema of all quantities stored
==========   ============================   ===============================================

.. _endianness: https://wizardzines.com/comics/little-endian/

The element by element data is stored in:

==========   ============================   ===============================================
``elem``     ``nel``-array of :ref:`elem`   array containing element data
==========   ============================   ===============================================


.. _elem:

Elem
----

This class contains the data for one hexahedral element organised as follows:

=========   ================   ====================================================
``pos``     array of floats    x,y,z coordinates of each grid point
``curv``    array of floats    radius of curvature of the element edges (*Nek5000*)
``ccurv``   array of strings   defines the type of curvature
``vel``     array of floats    u,v,w velocity at each grid point
``pres``    array of floats    pressure at each grid point
``temp``    array of floats    temperature at each grid point
``scal``    array of floats    passive scalars at each grid point
``bcs``     array of floats    list of boundary condition parameters (*Nek5000*)
=========   ================   ====================================================


.. _datalims:

DataLims
--------

This class contains the extrema of all quantities stored in the mesh.

========   ===============   ====================================================
``pos``    array of floats   max and min of x,y,z coordinates
``vel``    array of floats   max and min of u,v,w velocity
``pres``   array of floats   max and min of pressure
``temp``   array of floats   max and min of temperature
``scal``   array of floats   max and min of passive scalars
========   ===============   ====================================================



-------------------------------------------------------------------------------

.. _core-contents:

Contents of core.py
-------------------

.. automodule:: pymech.core
