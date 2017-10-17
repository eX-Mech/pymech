.. _exadata:
exadata.py
==========

This Python module contains the classes for the data structures used by
:ref:`neksuite` and :ref:`simsonsuite`.
The :ref:`exadata-contents` are reported at the bottom of this page.

exadata
-------

The main data class is ``exadata`` which is a structure for general and element
by element information on data stored on exahedral element(s).
The general information is constituted by:

==========   ============================   ===============================================
``ndim``     integer                        number of spatial dimensions of the data
``nel``      integer                        number of elements included in ``exadata``
``ncurv``    integer                        number of curved sides (used only by *Nek5000*)
``var``      string                         variables included in ``exadata`` (e.g. 'XUPT')
``lr1``      ``ndim``-array of integers     number of datapoints per element
``time``     float                          simulation time of the file
``istep``    integer                        simulation time step of the file
``wdsz``     integer                        word size, i.e. double or single precision file
``endian``   string                         endianness of the file (little/big)
``lims``     :ref:`datalims`                extrema of all quantities stored
==========   ============================   ===============================================

The element by element data is stored in:

==========   ============================   ===============================================
``elem``     ``nel``-array of :ref:`elem`   array containing element data
==========   ============================   ===============================================



.. _datalims:
datalims
--------


.. _elem:
elem
----

-------------------------------------------------------------------------------

.. _exadata-contents:
Contents of exadata.py
-----------------------

**exadata**

.. literalinclude:: ../src/exadata.py
   :language: python
   :lines: 53-69
