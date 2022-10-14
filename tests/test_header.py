import pytest

from pymech.neksuite.field import Header, read_header


def test_header_var_nb_vars_uninitialized():
    with pytest.raises(ValueError):
        Header(4, [8, 8, 8], 2, 2, 100, 10000, 3, 1)


def test_header_time(test_data_dir):
    fname = f"{test_data_dir}/nek/channel3D_0.f00001"
    header = read_header(fname)
    assert header.time == 0.2
