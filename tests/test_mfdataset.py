import copy
import shutil
from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def multi_data_dir(test_data_dir, tmp_path_factory):
    import pymech as pm

    file = test_data_dir / "nek" / "channel3D_0.f00001"
    fld = pm.readnek(file)
    # Copy to tmp_path_factory
    tmp_path = tmp_path_factory.mktemp("nek")
    shutil.copy(file, tmp_path)

    for step in range(2, 4):
        # In memory shallow copy
        new_fld = copy.copy(fld)
        # Note fld.time == 0.2
        new_fld.time = fld.time * step
        pm.writenek(tmp_path / f"channel3D_0.f0000{step}", new_fld)

    return tmp_path


def test_open_mfdataset_glob(multi_data_dir):
    from pymech.dataset import open_mfdataset

    ds = open_mfdataset(f"{multi_data_dir}/channel3D_0.f*")
    assert tuple(ds.time) == (0.2, 0.4, 0.6)


def test_open_mfdataset_files(multi_data_dir):
    from pymech.dataset import open_mfdataset

    files = sorted(Path(multi_data_dir).glob("channel3D_0.f*"))
    assert len(files) == 3

    ds = open_mfdataset(files)
    assert tuple(ds.time) == (0.2, 0.4, 0.6)


def test_xr_open_mfdataset(multi_data_dir):
    import xarray as xr

    ds = xr.open_mfdataset(
        f"{multi_data_dir}/channel3D_0.f*", combine="nested", concat_dim="time"
    )
    assert tuple(ds.time) == (0.2, 0.4, 0.6)
