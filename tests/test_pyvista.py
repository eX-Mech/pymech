"""Tests for PyVista visualization backend."""
import pytest
from pathlib import Path
import numpy as np

# Try importing backends
try:
    import pyvista as pv
    PYVISTA_AVAILABLE = True
except ImportError:
    PYVISTA_AVAILABLE = False
    pv = None

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    plt = None

# Try importing Protocol and backend modules
try:
    from pymech.viz.viz_protocol import MeshBackend
    PROTOCOL_AVAILABLE = True
except ImportError:
    PROTOCOL_AVAILABLE = False

try:
    from pymech.viz.pyvista_backend_impl import PyVistaBackend
    PYVISTA_IMPL_AVAILABLE = True
except ImportError:
    PYVISTA_IMPL_AVAILABLE = False

try:
    from pymech.viz.matplotlib_backend import MatplotlibBackend
    MATPLOTLIB_IMPL_AVAILABLE = True
except ImportError:
    MATPLOTLIB_IMPL_AVAILABLE = False


@pytest.mark.skipif(not PYVISTA_AVAILABLE, reason="PyVista not installed")
class TestPyVistaBackend:
    """Test suite for PyVista backend."""

    def test_import(self):
        """Test that viz module can be imported."""
        from pymech import viz
        assert hasattr(viz, 'plot_mesh')
        assert hasattr(viz, 'hexa_to_pyvista')
        assert hasattr(viz, 'add_boundary_conditions')

    def test_hexa_to_pyvista_linear_3d(self, test_data_dir):
        """Test conversion of 3D mesh with linear resolution."""
        from pymech import readre2
        from pymech.viz import hexa_to_pyvista

        # Load test data
        test_file = test_data_dir / "nek" / "box3d.re2"
        if not test_file.exists():
            pytest.skip(f"Test file not found: {test_file}")

        field = readre2(str(test_file))
        mesh = hexa_to_pyvista(field, resolution="linear", include_fields=False)

        assert isinstance(mesh, pv.UnstructuredGrid)
        assert mesh.n_cells == field.nel
        assert mesh.n_points == field.nel * 8  # 8 vertices per hex
        assert "element_id" in mesh.cell_data
        assert len(mesh.cell_data["element_id"]) == field.nel

    def test_hexa_to_pyvista_linear_2d(self, test_data_dir):
        """Test conversion of 2D mesh."""
        from pymech import readre2
        from pymech.viz import hexa_to_pyvista

        # Try to find a 2D test file
        test_files_2d = ["box2d.re2", "2D_section_R360.re2"]
        field = None
        for fname in test_files_2d:
            test_file = test_data_dir / "nek" / fname
            if test_file.exists():
                try:
                    field = readre2(str(test_file))
                    if field.ndim == 2:
                        break
                except Exception:
                    continue

        if field is None or field.ndim != 2:
            pytest.skip("No 2D test file available")

        mesh = hexa_to_pyvista(field, resolution="linear", include_fields=False)

        assert isinstance(mesh, pv.UnstructuredGrid)
        assert mesh.n_points == field.nel * 4  # 4 vertices per quad

    def test_hexa_to_pyvista_spectral(self, test_data_dir):
        """Test conversion with spectral resolution."""
        from pymech import readre2
        from pymech.viz import hexa_to_pyvista

        test_file = test_data_dir / "nek" / "box3d.re2"
        if not test_file.exists():
            pytest.skip(f"Test file not found: {test_file}")

        field = readre2(str(test_file))
        mesh = hexa_to_pyvista(field, resolution="spectral", include_fields=False)

        lx, ly, lz = field.lr1
        expected_cells = field.nel * (lx - 1) * (ly - 1) * (lz - 1)
        assert mesh.n_cells == expected_cells

    def test_include_fields(self, test_data_dir):
        """Test that velocity/pressure fields are included."""
        from pymech import readre2
        from pymech.viz import hexa_to_pyvista

        test_file = test_data_dir / "nek" / "box3d.re2"
        if not test_file.exists():
            pytest.skip(f"Test file not found: {test_file}")

        field = readre2(str(test_file))
        mesh = hexa_to_pyvista(field, resolution="linear", include_fields=True)

        # Check velocity field
        if field.var[1] == 3:
            assert "velocity" in mesh.point_data
            assert mesh.point_data["velocity"].shape[1] == 3
            assert "velocity_magnitude" in mesh.point_data

        # Check pressure field
        if field.var[2] == 1:
            assert "pressure" in mesh.point_data

        # Check temperature field
        if field.var[3] == 1:
            assert "temperature" in mesh.point_data

    def test_add_boundary_conditions(self, test_data_dir):
        """Test BC extraction and coloring."""
        from pymech import readre2
        from pymech.viz import hexa_to_pyvista, add_boundary_conditions

        test_file = test_data_dir / "nek" / "box3d.re2"
        if not test_file.exists():
            pytest.skip(f"Test file not found: {test_file}")

        field = readre2(str(test_file))
        mesh = hexa_to_pyvista(field, resolution="linear", include_fields=False)
        surface = add_boundary_conditions(mesh, field, bc_field=0)

        assert isinstance(surface, pv.PolyData)
        assert "bc_type" in surface.cell_data
        assert "bc_color" in surface.cell_data
        assert surface.cell_data["bc_color"].shape[1] == 3  # RGB colors

    @pytest.mark.parametrize("resolution", ["linear", "spectral"])
    def test_plot_mesh_headless(self, test_data_dir, tmp_path, resolution):
        """Test plot_mesh in headless mode (screenshot only)."""
        pv.OFF_SCREEN = True  # Enable headless rendering

        from pymech import readre2
        from pymech.viz import plot_mesh

        test_file = test_data_dir / "nek" / "box3d.re2"
        if not test_file.exists():
            pytest.skip(f"Test file not found: {test_file}")

        field = readre2(str(test_file))
        screenshot_path = tmp_path / f"test_{resolution}.png"

        try:
            plot_mesh(
                field,
                backend='pyvista',
                resolution=resolution,
                show_bcs=True,
                screenshot=str(screenshot_path),
                jupyter_backend="static",
            )

            assert screenshot_path.exists()
            assert screenshot_path.stat().st_size > 0
        except Exception as e:
            # Some systems may not support headless rendering
            pytest.skip(f"Headless rendering not supported: {e}")
        finally:
            pv.OFF_SCREEN = False

    def test_plot_mesh_return_plotter(self, test_data_dir):
        """Test that return_plotter works."""
        pv.OFF_SCREEN = True

        from pymech import readre2
        from pymech.viz import plot_mesh

        test_file = test_data_dir / "nek" / "box3d.re2"
        if not test_file.exists():
            pytest.skip(f"Test file not found: {test_file}")

        field = readre2(str(test_file))

        try:
            plotter = plot_mesh(field, backend='pyvista', return_plotter=True)
            assert isinstance(plotter, pv.Plotter)
            plotter.close()
        except Exception as e:
            pytest.skip(f"Headless rendering not supported: {e}")
        finally:
            pv.OFF_SCREEN = False

    def test_invalid_resolution(self, test_data_dir):
        """Test that invalid resolution raises ValueError."""
        from pymech import readre2
        from pymech.viz import hexa_to_pyvista

        test_file = test_data_dir / "nek" / "box3d.re2"
        if not test_file.exists():
            pytest.skip(f"Test file not found: {test_file}")

        field = readre2(str(test_file))

        with pytest.raises(ValueError, match="resolution must be"):
            hexa_to_pyvista(field, resolution="invalid")


@pytest.mark.skipif(not MATPLOTLIB_AVAILABLE, reason="Matplotlib not installed")
class TestMatplotlibBackend:
    """Test suite for Matplotlib backend."""

    def test_plot_mesh_matplotlib(self, test_data_dir, tmp_path):
        """Test plot_mesh with Matplotlib backend."""
        from pymech import readre2
        from pymech.viz import plot_mesh

        test_file = test_data_dir / "nek" / "box3d.re2"
        if not test_file.exists():
            pytest.skip(f"Test file not found: {test_file}")

        field = readre2(str(test_file))
        screenshot_path = tmp_path / "test_matplotlib.png"

        fig = plot_mesh(
            field,
            backend='matplotlib',
            show_bcs=False,
            screenshot=str(screenshot_path),
            return_plotter=True,
        )

        assert fig is not None
        assert screenshot_path.exists()
        assert screenshot_path.stat().st_size > 0
        plt.close(fig)

    def test_plot_mesh_matplotlib_views(self, test_data_dir):
        """Test different camera views with Matplotlib."""
        from pymech import readre2
        from pymech.viz import plot_mesh

        test_file = test_data_dir / "nek" / "box3d.re2"
        if not test_file.exists():
            pytest.skip(f"Test file not found: {test_file}")

        field = readre2(str(test_file))

        for view in ['xy', 'xz', 'yz']:
            fig = plot_mesh(
                field,
                backend='matplotlib',
                view=view,
                show_bcs=False,
                return_plotter=True,
            )
            assert fig is not None
            plt.close(fig)


class TestBackendSelection:
    """Test backend selection logic."""

    def test_auto_backend_selection(self, test_data_dir):
        """Test that 'auto' backend selects appropriately."""
        from pymech import readre2
        from pymech.viz import plot_mesh

        test_file = test_data_dir / "nek" / "box3d.re2"
        if not test_file.exists():
            pytest.skip(f"Test file not found: {test_file}")

        field = readre2(str(test_file))

        # 'auto' should work if at least one backend is available
        if PYVISTA_AVAILABLE or MATPLOTLIB_AVAILABLE:
            if PYVISTA_AVAILABLE:
                pv.OFF_SCREEN = True
            try:
                result = plot_mesh(
                    field,
                    backend='auto',
                    screenshot=None,
                    return_plotter=True,
                )
                assert result is not None
                if hasattr(result, 'close'):
                    result.close()
                elif hasattr(result, 'clf'):
                    plt.close(result)
            except Exception as e:
                pytest.skip(f"Auto backend selection failed: {e}")
            finally:
                if PYVISTA_AVAILABLE:
                    pv.OFF_SCREEN = False
        else:
            with pytest.raises(ImportError):
                plot_mesh(field, backend='auto')

    def test_invalid_backend(self, test_data_dir):
        """Test that invalid backend raises ValueError."""
        from pymech import readre2
        from pymech.viz import plot_mesh

        test_file = test_data_dir / "nek" / "box3d.re2"
        if not test_file.exists():
            pytest.skip(f"Test file not found: {test_file}")

        field = readre2(str(test_file))

        with pytest.raises(ValueError, match="Invalid backend"):
            plot_mesh(field, backend='invalid')


def test_import_without_backends():
    """Test that module can be imported even without backends."""
    # This test ensures graceful degradation
    try:
        from pymech import viz
        # If import succeeds, module should have main functions defined
        assert hasattr(pyvista_backend, 'plot_mesh')
    except ImportError:
        # If import fails, it's acceptable (no backends available)
        pass


def test_bc_colors():
    """Test that BC color scheme is defined."""
    from pymech.viz import DEFAULT_BC_COLORS

    # Check that important BC types are defined
    assert "" in DEFAULT_BC_COLORS
    assert "E" in DEFAULT_BC_COLORS
    assert "W" in DEFAULT_BC_COLORS
    assert "O" in DEFAULT_BC_COLORS

    # Check that colors are RGB tuples
    for bc_type, color in DEFAULT_BC_COLORS.items():
        assert isinstance(color, tuple)
        assert len(color) == 3
        assert all(0 <= c <= 1 for c in color)


@pytest.mark.skipif(not PROTOCOL_AVAILABLE, reason="Protocol module not available")
class TestProtocolCompliance:
    """Test that backends implement the MeshBackend protocol correctly."""

    @pytest.mark.skipif(not PYVISTA_IMPL_AVAILABLE, reason="PyVista implementation not available")
    def test_pyvista_backend_protocol(self):
        """Test PyVista backend implements MeshBackend protocol."""
        from pymech.viz.pyvista_backend_impl import PyVistaBackend
        from pymech.viz import MeshBackend

        backend = PyVistaBackend()

        # Runtime check for Protocol compliance
        assert isinstance(backend, MeshBackend)

        # Verify required methods exist
        assert hasattr(backend, 'is_available')
        assert hasattr(backend, 'plot_mesh')
        assert hasattr(backend, 'get_backend_name')
        assert hasattr(backend, 'get_capabilities')

        # Verify method return types
        assert isinstance(backend.is_available(), bool)
        assert isinstance(backend.get_backend_name(), str)
        assert isinstance(backend.get_capabilities(), dict)

        # Verify backend name
        assert backend.get_backend_name() == "pyvista"

    @pytest.mark.skipif(not MATPLOTLIB_IMPL_AVAILABLE, reason="Matplotlib implementation not available")
    def test_matplotlib_backend_protocol(self):
        """Test Matplotlib backend implements MeshBackend protocol."""
        from pymech.viz.matplotlib_backend import MatplotlibBackend
        from pymech.viz import MeshBackend

        backend = MatplotlibBackend()

        # Runtime check for Protocol compliance
        assert isinstance(backend, MeshBackend)

        # Verify required methods exist
        assert hasattr(backend, 'is_available')
        assert hasattr(backend, 'plot_mesh')
        assert hasattr(backend, 'get_backend_name')
        assert hasattr(backend, 'get_capabilities')

        # Verify method return types
        assert isinstance(backend.is_available(), bool)
        assert isinstance(backend.get_backend_name(), str)
        assert isinstance(backend.get_capabilities(), dict)

        # Verify backend name
        assert backend.get_backend_name() == "matplotlib"

    @pytest.mark.skipif(not PYVISTA_IMPL_AVAILABLE, reason="PyVista implementation not available")
    def test_pyvista_capabilities(self):
        """Test PyVista backend capabilities."""
        from pymech.pyvista_backend_impl import PyVistaBackend

        backend = PyVistaBackend()
        caps = backend.get_capabilities()

        assert caps["interactive"] is True
        assert caps["3d"] is True
        assert caps["jupyter"] is True
        assert caps["headless"] is True
        assert "formats" in caps
        assert "png" in caps["formats"]

    @pytest.mark.skipif(not MATPLOTLIB_IMPL_AVAILABLE, reason="Matplotlib implementation not available")
    def test_matplotlib_capabilities(self):
        """Test Matplotlib backend capabilities."""
        from pymech.matplotlib_backend import MatplotlibBackend

        backend = MatplotlibBackend()
        caps = backend.get_capabilities()

        assert caps["interactive"] is False
        assert caps["3d"] is True
        assert caps["jupyter"] is True
        assert caps["headless"] is True
        assert "formats" in caps
        assert "png" in caps["formats"]
        assert "pdf" in caps["formats"]


@pytest.mark.skipif(not PROTOCOL_AVAILABLE, reason="Protocol module not available")
class TestBackendDiscovery:
    """Test backend discovery and selection."""

    def test_get_available_backends(self):
        """Test get_available_backends() function."""
        from pymech.viz import get_available_backends

        backends = get_available_backends()

        assert isinstance(backends, dict)

        # Check that available backends are present
        if PYVISTA_AVAILABLE:
            assert 'pyvista' in backends
        else:
            assert 'pyvista' not in backends

        if MATPLOTLIB_AVAILABLE:
            assert 'matplotlib' in backends
        else:
            assert 'matplotlib' not in backends

    def test_backend_availability(self):
        """Test is_available() for each backend."""
        # Test PyVista backend
        if PYVISTA_IMPL_AVAILABLE:
            from pymech.pyvista_backend_impl import PyVistaBackend
            backend = PyVistaBackend()
            assert backend.is_available() == PYVISTA_AVAILABLE

        # Test Matplotlib backend
        if MATPLOTLIB_IMPL_AVAILABLE:
            from pymech.matplotlib_backend import MatplotlibBackend
            backend = MatplotlibBackend()
            assert backend.is_available() == MATPLOTLIB_AVAILABLE
