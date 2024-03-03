"""Test xarray dataset interface"""

import math

import numpy.testing as npt
import pytest
import xarray as xr


def _check_channel_dset(ds):
    assert isinstance(ds, xr.Dataset)
    assert all(
        var in ds.data_vars for var in ("ux", "uy", "uz", "xmesh", "ymesh", "zmesh")
    )
    assert all(coord in ds.coords for coord in ("x", "y", "z"))
    assert tuple(ds.dims.values()) == (64, 64, 64)

    assert math.isclose(ds.x.max(), 2 * math.pi, abs_tol=1e-6)
    assert math.isclose(ds.y.max(), 1.0, abs_tol=1e-6)
    assert math.isclose(ds.z.max(), math.pi, abs_tol=1e-6)
    assert math.isclose(ds.time, 0.2, abs_tol=1e-6)


def test_channel_dataset(test_data_dir):
    import pymech.dataset as pd

    fname = f"{test_data_dir}/nek/channel3D_0.f00001"
    ds = pd.open_dataset(fname)
    _check_channel_dset(ds)


def test_phill_dataset(test_data_dir):
    import pymech.dataset as pd

    fname = f"{test_data_dir}/nek/phill0.f00010"
    with pytest.raises(
        NotImplementedError, match="currently works only with cartesian box meshes"
    ):
        pd.open_dataset(fname)


@pytest.mark.parametrize("file_name", ["channel3D_0.f00001", "stsabl0.f00001"])
def test_dataset_coords(test_data_dir, file_name):
    """Check if the 1D coordinates match with the original 3D coordinate arrays"""
    import pymech.dataset as pd

    path = f"{test_data_dir}/nek/" + file_name
    ds = pd.open_dataset(path)

    dsx = ds.mean(("y", "z"))
    dsy = ds.mean(("x", "z"))

    npt.assert_allclose(ds.x.data, dsx.xmesh.data, rtol=1e-6)
    npt.assert_allclose(ds.y.data, dsy.ymesh.data, rtol=1e-6)


def test_nek_ext_regex():
    from pymech.dataset import nek_ext_pattern as pattern

    assert not pattern.match(".f90")
    assert not pattern.match(".f")
    assert not pattern.match(".fort")
    assert not pattern.match(".f0000")
    assert not pattern.match(".fll")
    assert pattern.match(".fld")
    assert pattern.match(".f00001")
    assert pattern.match("case0.f12345")


def test_xarray_plugin_entrypoint():
    plugins = xr.backends.plugins
    assert "pymech" in plugins.list_engines()
    assert plugins.get_backend("pymech")


def test_xarray_open_dataset(test_data_dir):
    # With engine explicitly mentioned
    xr.open_dataset(f"{test_data_dir}/nek/channel3D_0.f00001", engine="pymech")

    # Let xarray guess the engine
    ds = xr.open_dataset(f"{test_data_dir}/nek/channel3D_0.f00001")
    _check_channel_dset(ds)
