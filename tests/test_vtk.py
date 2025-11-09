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
@pytest.mark.xfail(
    raises=ValueError, reason="VTK is built using numpy 1.x and produces ValueError"
)
@pytest.mark.parametrize("downsample", (True, False))
def test_writevtk(downsample, test_data_dir, tmpdir):
    from pymech import readnek
    from pymech.vtksuite import writevtk

    in_fname = Path(test_data_dir) / "nek" / "channel3D_0.f00001"
    out_fname = Path(tmpdir / in_fname.name)
    field = readnek(in_fname)

    import traits

    try:
        writevtk(out_fname, field)
    except traits.trait_errors.TraitError:
        pytest.xfail(
            reason="Known issue with changes in tvtk. See https://github.com/eX-Mech/pymech/issues/142"
        )
    assert out_fname.with_suffix(".vtp").exists()
