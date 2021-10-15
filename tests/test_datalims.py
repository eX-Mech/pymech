import pytest
from numpy import testing as npt


def test_lims(test_data_dir):
    from pymech.neksuite import readnek

    fld = readnek(f"{test_data_dir}/nek/channel3D_0.f00001")

    with pytest.raises(AttributeError):
        fld.lims = ()

    with pytest.raises(AttributeError):
        fld.lims.pos = ((0, 3.14), (0, 1), (-1, 1))

    with pytest.raises(TypeError):
        fld.lims.pos[0, 1] = 3.14

    npt.assert_array_equal(
        fld.lims.pos,
        (
            (0.0, 6.2831854820251465),
            (-1.0, 1.0),
            (0.0, 3.1415927410125732),
        ),
    )
    npt.assert_array_equal(
        fld.lims.vel,
        (
            (-0.0001499584031989798, 1.2596282958984375),
            (-0.024550313130021095, 0.021655898541212082),
            (-0.022412149235606194, 0.020889850333333015),
        ),
    )
    npt.assert_array_equal(fld.lims.pres, ((-0.2145293802022934, 0.2344336062669754),))
    npt.assert_array_equal(fld.lims.temp, ())
    npt.assert_array_equal(fld.lims.scal, ())
