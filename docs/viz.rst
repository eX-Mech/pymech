.. _viz:

pymech.viz
==========

3D mesh visualization using PyVista or Matplotlib backends.

.. warning::

   This subpackage requires visualization libraries. Install with::

       pip install pymech[plot]

   This installs PyVista, Matplotlib, and required dependencies.

Overview
--------

The ``pymech.viz`` subpackage provides modern mesh visualization capabilities with two backends:

- **PyVista** (recommended): Interactive 3D visualization optimized for Jupyter notebooks
- **Matplotlib**: Publication-quality figures and simple 3D plots

Architecture
~~~~~~~~~~~~

The visualization system uses a modular, protocol-based architecture within the ``pymech.viz`` subpackage:

- ``pymech.viz.viz_protocol``: Defines the ``MeshBackend`` Protocol interface
- ``pymech.viz.pyvista_backend_impl``: PyVista-specific implementation
- ``pymech.viz.matplotlib_backend``: Matplotlib-specific implementation
- ``pymech.viz.pyvista_backend``: Backend dispatcher module

The subpackage exports a unified API through ``pymech.viz``:

- Main functions: ``plot_mesh()``, ``get_available_backends()``
- PyVista utilities: ``hexa_to_pyvista()``, ``add_boundary_conditions()``
- Protocol and constants: ``MeshBackend``, ``DEFAULT_BC_COLORS``

Backends implement the ``MeshBackend`` protocol using ``typing.Protocol`` and
``runtime_checkable``, ensuring consistent interfaces across different rendering engines.
The main ``plot_mesh()`` function automatically selects the best available backend or
uses the explicitly specified one

Quick Start
-----------

PyVista Backend (Interactive 3D)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import pymech as pm
    from pymech.viz import plot_mesh

    # Load mesh
    field = pm.readnek("channel3D_0.f00001")

    # Plot with boundary conditions
    plot_mesh(field, backend='pyvista')

Matplotlib Backend (Publication Figures)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    # Plot with Matplotlib for static figure
    plot_mesh(field, backend='matplotlib', view='xy',
              screenshot='mesh_figure.png')

Jupyter Notebook Usage
----------------------

PyVista automatically detects Jupyter environments and uses interactive backends:

.. code-block:: python

    # In Jupyter notebook - automatically interactive
    from pymech.viz import plot_mesh

    field = pm.readnek("mesh.nek5000")
    plot_mesh(field, jupyter_backend="trame")

Available Jupyter backends:

- ``"trame"``: Interactive (default, requires trame package)
- ``"static"``: Non-interactive PNG/JPEG
- ``"ipyvtklink"``: Alternative interactive backend

Features
--------

Mesh Resolution Modes
~~~~~~~~~~~~~~~~~~~~~

**Linear resolution** (default, fast):

.. code-block:: python

    plot_mesh(field, resolution='linear')

Uses only corner vertices (8 per hexahedron), suitable for most visualizations.

**Spectral resolution** (accurate, slower):

.. code-block:: python

    plot_mesh(field, resolution='spectral')

Uses all GLL points, subdividing each element. Preserves curved geometries accurately.

Boundary Condition Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Color-coded boundary conditions:

.. code-block:: python

    # Show BCs (default)
    plot_mesh(field, show_bcs=True, bc_field=0)

    # bc_field=0: velocity BCs
    # bc_field=1: temperature BCs

BC Color Scheme:

- **W** (wall): Dark blue
- **v** (velocity BC): Light blue
- **O** (outflow): Dark red
- **T** (temperature): Dark green
- **E** (element): Black
- **P** (periodic): Gray

Field Visualization
~~~~~~~~~~~~~~~~~~~

Visualize velocity, pressure, or temperature fields:

.. code-block:: python

    from pymech.viz import hexa_to_pyvista
    import pyvista as pv

    # Convert mesh with field data
    mesh = hexa_to_pyvista(field, include_fields=True)

    # Plot velocity magnitude
    plotter = pv.Plotter()
    plotter.add_mesh(mesh, scalars="velocity_magnitude",
                     cmap="coolwarm", show_edges=True)
    plotter.show()

Custom Visualization
--------------------

Advanced customization with PyVista:

.. code-block:: python

    # Return plotter for customization
    plotter = plot_mesh(field, return_plotter=True)

    # Customize camera
    plotter.camera_position = [(10, 10, 10), (0, 0, 0), (0, 1, 0)]

    # Add text annotation
    plotter.add_text("Channel Flow Mesh", position="upper_right",
                     font_size=12)

    # Show
    plotter.show()

Export Screenshots
~~~~~~~~~~~~~~~~~~

.. code-block:: python

    # Export high-resolution image
    plot_mesh(field, screenshot="mesh_hires.png")

    # Matplotlib for publication figures
    plot_mesh(field, backend='matplotlib',
              screenshot="mesh_pub.pdf", figsize=(12, 10))

Backend Comparison
------------------

+----------------------+---------------------+----------------------+
| Feature              | PyVista             | Matplotlib           |
+======================+=====================+======================+
| Dimensionality       | Full 3D             | Limited 3D           |
+----------------------+---------------------+----------------------+
| Interactivity        | Full rotation/zoom  | Basic                |
+----------------------+---------------------+----------------------+
| Jupyter support      | Excellent (trame)   | Good                 |
+----------------------+---------------------+----------------------+
| BC visualization     | Face colors         | Edge colors          |
+----------------------+---------------------+----------------------+
| Performance          | Good for large mesh | Slow for large mesh  |
+----------------------+---------------------+----------------------+
| Use case             | Exploration         | Publication figures  |
+----------------------+---------------------+----------------------+

**Recommendation**: Use PyVista for interactive exploration and Jupyter workflows.
Use Matplotlib for generating publication-quality 2D projections.

API Reference
-------------

Main Function
~~~~~~~~~~~~~

.. autofunction:: pymech.viz.plot_mesh
   :noindex:

Backend Discovery
~~~~~~~~~~~~~~~~~

.. autofunction:: pymech.viz.get_available_backends
   :noindex:

Conversion Functions (PyVista)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: pymech.viz.hexa_to_pyvista
   :noindex:

.. autofunction:: pymech.viz.add_boundary_conditions
   :noindex:

Complete API
~~~~~~~~~~~~

.. automodule:: pymech.viz
   :members:
   :undoc-members:
   :show-inheritance:

Examples Gallery
----------------

Example 1: Check Available Backends
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from pymech.viz import get_available_backends

    backends = get_available_backends()
    print(f"Available backends: {list(backends.keys())}")

    # Check capabilities
    for name, backend in backends.items():
        caps = backend.get_capabilities()
        print(f"{name}: 3D={caps['3d']}, Interactive={caps['interactive']}")

Example 2: Basic Mesh Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import pymech as pm
    from pymech.viz import plot_mesh

    field = pm.readnek("channel3D_0.f00001")
    plot_mesh(field)

Example 3: Custom Camera View
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    plotter = plot_mesh(field, return_plotter=True)
    plotter.camera_position = 'xy'  # Top-down view
    plotter.show()

Example 4: Velocity Field Coloring
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from pymech.viz import hexa_to_pyvista
    import pyvista as pv

    mesh = hexa_to_pyvista(field, include_fields=True)
    mesh.plot(scalars="pressure", cmap="viridis", show_edges=True)

Example 5: Publication Figure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    plot_mesh(field, backend='matplotlib', view='xz',
              show_bcs=False, color='darkblue',
              screenshot='figure.pdf', figsize=(8, 6))

Comparison with meshplot
-------------------------

The existing ``meshplot.py`` module provides 2D visualization with wxPython.
Choose the appropriate tool for your use case:

+-----------------------+---------------------+---------------------------+
| Feature               | meshplot (2D)       | pymech.viz (3D)           |
+=======================+=====================+===========================+
| Dimensionality        | 2D only             | 2D and 3D                 |
+-----------------------+---------------------+---------------------------+
| Curved edges          | Exact (parabolic)   | Linear approx.            |
+-----------------------+---------------------+---------------------------+
| Jupyter support       | No                  | Yes (primary)             |
+-----------------------+---------------------+---------------------------+
| Boundary conditions   | Edge colors         | Face colors               |
+-----------------------+---------------------+---------------------------+
| Interactivity         | Zoom/pan only       | Full 3D rotation          |
+-----------------------+---------------------+---------------------------+
| Dependencies          | wxPython, OpenGL    | PyVista or Matplotlib     |
+-----------------------+---------------------+---------------------------+

**Use ``meshplot``** for detailed 2D edge inspection with exact curved geometry.

**Use ``pymech.viz``** for 3D exploration, Jupyter workflows, and publication figures.

Troubleshooting
---------------

PyVista Not Available
~~~~~~~~~~~~~~~~~~~~~

If you get "PyVista is required" error:

.. code-block:: bash

    pip install pymech[plot]

or install PyVista separately:

.. code-block:: bash

    pip install pyvista trame ipywidgets

Headless Rendering
~~~~~~~~~~~~~~~~~~

For headless servers or CI/CD:

.. code-block:: python

    import pyvista as pv
    pv.OFF_SCREEN = True

    plot_mesh(field, screenshot="output.png")

Jupyter Display Issues
~~~~~~~~~~~~~~~~~~~~~~

If plots don't show in Jupyter:

.. code-block:: python

    import pyvista as pv
    pv.set_jupyter_backend('static')  # Try static backend

    # Or try different backend
    plot_mesh(field, jupyter_backend='ipyvtklink')

See Also
--------

- :ref:`vtksuite` - VTK export functions
- :ref:`meshtools` - Mesh manipulation utilities
- :ref:`dataset` - Xarray integration
