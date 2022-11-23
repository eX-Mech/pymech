"""Tests experimental vtksuite"""
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
def test_writevtk(test_data_dir, downsample):
    from pymech import readnek
    from pymech.vtksuite import writevtk

    fname = f"{test_data_dir}/nek/channel3D_0.f00001"
    field = readnek(fname)
    writevtk(fname, field)
    assert fname.with_suffix(".vtp").exists()
