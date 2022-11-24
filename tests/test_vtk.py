"""Tests experimental vtksuite"""
from pathlib import Path
import pytest

try:
    import tvtk  # noqa
except ImportError:
    TVTK_INSTALLED = False
else:
    TVTK_INSTALLED = True


@pytest.mark.skipif(
    not TVTK_INSTALLED, reason="Package mayavi / tvtk is not installed "
)
@pytest.mark.parametrize("downsample", (True, False))
def test_writevtk(downsample, test_data_dir, tmpdir):
    from pymech import readnek
    from pymech.vtksuite import writevtk

    in_fname = Path(test_data_dir) / "nek" / "channel3D_0.f00001"
    out_fname = Path(tmpdir / in_fname.name)
    field = readnek(in_fname)
    writevtk(out_fname, field)
    assert out_fname.with_suffix(".vtp").exists()
