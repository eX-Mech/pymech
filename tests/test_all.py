import os
import logging
import math
import time
import tempfile
from pathlib import Path
from textwrap import dedent

from numpy import testing as npt
from pymech.log import logger


def setup_module(module):
    logger.setLevel(logging.DEBUG)

    tmp = Path(tempfile.mkdtemp(suffix="pymech_tests"))
    (tmp / "tests").symlink_to(Path(__file__).parent)
    os.chdir(tmp)

    logger.info(f"Test directory: {tmp}")


# ------------------------------------------------------------------------------
# test nek scripts
#
def test_readnek():
    import pymech.neksuite as ns

    fname = "./tests/nek/channel3D_0.f00001"
    field = ns.readnek(fname)

    assert field.endian == "little"
    assert field.istep == 10
    assert field.lr1 == [8, 8, 8]
    assert field.ndim == 3
    assert field.nel == 512
    assert field.var == [3, 3, 1, 0, 0]
    assert field.wdsz == 4
    assert (field.time - 0.2) < 1e-3
    representation = dedent(
        """\
    <pymech.exadata.exadata>
    Dimensions:    3
    Precision:     4 bytes
    Mesh limits:
      * x:         [0.         6.28318548]
      * y:         [-1.  1.]
      * z:         [0.         3.14159274]
    Time:
      * time:      0.2
      * istep:     10
    Elements:
      * nel:       512
      * elem:      [<elem centered at [ 0.39269908 -0.98        0.19634954]>
                    ...
                    <elem centered at [5.89048618 0.98       2.94524309]>]"""
    )
    assert repr(field).splitlines() == representation.splitlines()


def test_writenek():
    import pymech.neksuite as ns

    fname = "./tests/nek/channel3D_0.f00001"
    time0 = time.perf_counter()
    field = ns.readnek(fname)
    time1 = time.perf_counter()

    fnamew = "./test_0.f00001"
    status = ns.writenek(fnamew, field)
    time2 = time.perf_counter()
    logger.info(
        "readnek: {:.6e} s; writenek: {:.6e} s".format(time1 - time0, time2 - time1)
    )

    assert status == 0

    fieldw = ns.readnek(fnamew)

    assert field.endian == fieldw.endian
    assert field.istep == fieldw.istep
    assert field.lr1 == fieldw.lr1
    assert field.ndim == fieldw.ndim
    assert field.nel == fieldw.nel
    assert field.var == fieldw.var
    assert field.wdsz == fieldw.wdsz
    assert (field.time - fieldw.time) < 1e-3
    npt.assert_array_equal(field.lims.pos, fieldw.lims.pos)
    npt.assert_array_equal(field.lims.vel, fieldw.lims.vel)
    npt.assert_array_equal(field.lims.pres, fieldw.lims.pres)
    npt.assert_array_equal(field.lims.scal, fieldw.lims.scal)

    for elem, elemw in zip(field.elem, fieldw.elem):
        npt.assert_array_equal(elem.pos, elemw.pos)
        npt.assert_array_equal(elem.vel, elemw.vel)
        npt.assert_array_equal(elem.pres, elemw.pres)
        npt.assert_array_equal(elem.scal, elemw.scal)


def test_writenek_big_endian():
    import pymech.neksuite as ns

    fname = "./tests/nek/channel3D_0.f00001"
    field = ns.readnek(fname)

    field.endian = "big"
    fnamew = "./test_0_big.f00001"
    status = ns.writenek(fnamew, field)
    assert status == 0

    fieldw = ns.readnek(fnamew)
    assert field.endian == fieldw.endian == "big"

    npt.assert_array_equal(field.lims.pos, fieldw.lims.pos)
    npt.assert_array_equal(field.lims.vel, fieldw.lims.vel)
    npt.assert_array_equal(field.lims.pres, fieldw.lims.pres)
    npt.assert_array_equal(field.lims.scal, fieldw.lims.scal)

    for elem, elemw in zip(field.elem, fieldw.elem):
        npt.assert_array_equal(elem.pos, elemw.pos)
        npt.assert_array_equal(elem.vel, elemw.vel)
        npt.assert_array_equal(elem.pres, elemw.pres)
        npt.assert_array_equal(elem.scal, elemw.scal)


def test_readnek_scalars():
    import pymech.neksuite as ns

    # 2D statistics file
    fname = "./tests/nek/stsabl0.f00001"
    field = ns.readnek(fname)

    ux_min, ux_max = field.lims.scal[0]
    assert math.isclose(ux_max, 5.3, abs_tol=0.1)


def test_writenek_scalars():
    import pymech.neksuite as ns

    fname = "./tests/nek/stsabl0.f00001"
    field = ns.readnek(fname)

    fnamew = "./test_sts_0.f00001"
    status = ns.writenek(fnamew, field)

    assert status == 0

    fieldw = ns.readnek(fnamew)
    npt.assert_array_equal(field.lims.scal, fieldw.lims.scal)


def test_readrea():
    import pymech.neksuite as ns

    fname = "./tests/nek/2D_section_R360.rea"
    field = ns.readrea(fname)

    assert field.lr1 == [2, 2, 1]
    assert field.ndim == 2
    assert field.nel == 1248
    assert abs(field.elem[0].pos[0][0][0][0] - 0.048383219999999998) < 1e-3
    assert abs(field.elem[887].curv[1, 0] - 1.21664) < 1e-3
    assert field.elem[887].ccurv[1] == "C"

    fname = "./tests/nek/m3j_bf_test.rea"
    field = ns.readrea(fname)
    assert field.elem[790].ccurv[0] == "m"
    assert abs(field.elem[790].curv[0][1] + 0.05258981) < 1e-7
    assert field.elem[0].bcs[0, 0][0] == "W"
    assert field.elem[0].bcs[0, 1][0] == "o"
    assert field.elem[0].bcs[0, 2][0] == "E"
    assert field.elem[0].bcs[0, 2][1] == 1
    assert field.elem[0].bcs[0, 2][2] == 3
    assert int(field.elem[0].bcs[0, 2][3]) == 2
    assert int(field.elem[0].bcs[0, 2][4]) == 1
    assert int(field.elem[799].bcs[1, 1][3]) == 790
    assert field.elem[799].bcs[1, 2][0] == "t"
    assert field.elem[799].bcs[1, 3][0] == "I"
    assert int(field.elem[799].bcs[2, 1][3]) == 790
    assert field.elem[799].bcs[2, 2][0] == "P"
    assert field.elem[799].bcs[2, 3][0] == "P"


def test_writerea():
    import pymech.neksuite as ns

    fname = "./tests/nek/2D_section_R360.rea"
    field = ns.readrea(fname)

    fnamew = "test.rea"
    status = ns.writerea(fnamew, field)

    assert status == 0

    fieldw = ns.readrea(fnamew)

    assert field.endian == fieldw.endian
    assert field.lr1 == fieldw.lr1
    assert field.ndim == fieldw.ndim
    assert field.nel == fieldw.nel
    assert field.wdsz == fieldw.wdsz
    assert (field.elem[0].pos[0][0][0][0] - fieldw.elem[0].pos[0][0][0][0]) < 1e-3
    assert abs(field.elem[887].curv[1, 0] - 1.21664) < 1e-3
    assert field.elem[887].ccurv[1] == "C"

    fname = "./tests/nek/m3j_bf_test.rea"
    fnamew = "test.rea"

    field = ns.readrea(fname)
    status = ns.writerea(fnamew, field)

    assert status == 0

    fieldw = ns.readrea(fnamew)

    assert fieldw.elem[790].ccurv[0] == "m"
    assert abs(fieldw.elem[790].curv[0][1] + 0.05258981) < 1e-7
    assert fieldw.elem[0].bcs[0, 0][0] == "W"
    assert fieldw.elem[0].bcs[0, 1][0] == "o"
    assert fieldw.elem[0].bcs[0, 2][0] == "E"
    assert fieldw.elem[0].bcs[0, 2][1] == 1
    assert fieldw.elem[0].bcs[0, 2][2] == 3
    assert int(fieldw.elem[0].bcs[0, 2][3]) == 2
    assert int(fieldw.elem[0].bcs[0, 2][4]) == 1
    assert int(field.elem[799].bcs[1, 1][3]) == 790
    assert fieldw.elem[799].bcs[1, 2][0] == "t"
    assert fieldw.elem[799].bcs[1, 3][0] == "I"
    assert int(fieldw.elem[799].bcs[2, 1][3]) == 790
    assert fieldw.elem[799].bcs[2, 2][0] == "P"
    assert fieldw.elem[799].bcs[2, 3][0] == "P"


def test_merge():
    import pymech.neksuite as ns
    import copy

    fname = "./tests/nek/box3d.rea"
    mesh = ns.readrea(fname)
    mesh1 = copy.deepcopy(mesh)
    mesh2 = copy.deepcopy(mesh)
    mesh3 = copy.deepcopy(mesh)
    # mesh and mesh1 will be connected along the 'v' and 'O' BCs
    for el in mesh1.elem:
        el.pos[0, ...] = el.pos[0, ...] + 2
    # mesh and mesh2 will be connected along 'P' BCs
    for el in mesh2.elem:
        el.pos[1, ...] = el.pos[1, ...] + 2
    for el in mesh3.elem:
        el.pos[2, ...] = el.pos[2, ...] + 2
    n1 = mesh1.merge(mesh)
    n2 = mesh2.merge(mesh)
    n3 = mesh3.merge(mesh)
    assert mesh.check_connectivity()
    assert mesh1.check_connectivity()
    assert mesh2.check_connectivity()
    assert mesh3.check_connectivity()
    assert mesh1.nel == 2 * mesh.nel
    assert n1 == 9
    assert n2 == 9
    assert n3 == 9
    # check if the element/faces indices in the boundary conditions are right too, even if it may not matter
    assert mesh1.nbc > 0
    for ibc in range(mesh1.nbc):
        for (iel, el) in enumerate(mesh1.elem):
            for iface in range(6):
                assert el.bcs[ibc, iface][1] == iel + 1
                assert el.bcs[ibc, iface][2] == iface + 1


def test_readre2():
    import pymech.neksuite as ns

    # The .re2 has been generated with reatore2 and contains the same data
    # except for the internal boundary conditions.
    # Assuming that `readrea` is correct, this checks id the .re2 file is read correctly too.
    #
    frea = "./tests/nek/2D_section_R360.rea"
    fre2 = "./tests/nek/2D_section_R360.re2"
    meshrea = ns.readrea(frea)
    meshre2 = ns.readre2(fre2)

    # remove the 'E' conditions from the .rea data
    for el in meshrea.elem:
        for iface in range(4):
            if el.bcs[0, iface][0] == "E":
                el.bcs[0, iface][0] = ""
                for j in range(1, 8):
                    el.bcs[0, iface][j] = 0

    assert meshre2.ndim == meshrea.ndim
    assert meshre2.nel == meshrea.nel
    assert meshre2.ncurv == meshrea.ncurv
    assert meshre2.nbc == meshrea.nbc
    assert meshre2.var == meshrea.var
    assert meshre2.lr1 == meshrea.lr1
    assert meshre2.wdsz == 8
    for (el, elw) in zip(meshrea.elem, meshre2.elem):
        npt.assert_allclose(elw.pos, el.pos)
        npt.assert_array_equal(elw.bcs, el.bcs)
        npt.assert_allclose(elw.curv, el.curv)
        npt.assert_array_equal(elw.ccurv, el.ccurv)


def test_readre2_3d():
    import pymech.neksuite as ns

    # same test as test_readre2(), but with a 3D mesh.

    frea = "./tests/nek/box3d.rea"
    fre2 = "./tests/nek/box3d.re2"
    meshrea = ns.readrea(frea)
    meshre2 = ns.readre2(fre2)
    # remove the 'E' conditions from the .rea data
    for el in meshrea.elem:
        for iface in range(6):
            if el.bcs[0, iface][0] == "E":
                el.bcs[0, iface][0] = ""
                for j in range(1, 8):
                    el.bcs[0, iface][j] = 0

    assert meshre2.ndim == meshrea.ndim
    assert meshre2.nel == meshrea.nel
    assert meshre2.ncurv == meshrea.ncurv
    assert meshre2.nbc == meshrea.nbc
    assert meshre2.var == meshrea.var
    assert meshre2.lr1 == meshrea.lr1
    assert meshre2.wdsz == 8
    for (el, elw) in zip(meshrea.elem, meshre2.elem):
        npt.assert_allclose(elw.pos, el.pos)
        npt.assert_array_equal(elw.bcs, el.bcs)
        npt.assert_allclose(elw.curv, el.curv)
        npt.assert_array_equal(elw.ccurv, el.ccurv)


def test_writere2():
    import pymech.neksuite as ns

    fin = "./tests/nek/2D_section_R360.re2"
    fout = "./test_1.re2"
    mesh = ns.readre2(fin)
    status = ns.writere2(fout, mesh)

    assert status == 0

    meshw = ns.readre2(fout)

    assert meshw.ndim == mesh.ndim
    assert meshw.nel == mesh.nel
    assert meshw.ncurv == mesh.ncurv
    assert meshw.nbc == mesh.nbc
    assert meshw.var == mesh.var
    assert meshw.lr1 == mesh.lr1
    assert meshw.wdsz == 8
    for (el, elw) in zip(mesh.elem, meshw.elem):
        npt.assert_array_equal(elw.pos, el.pos)
        npt.assert_array_equal(elw.bcs, el.bcs)
        npt.assert_array_equal(elw.curv, el.curv)
        npt.assert_array_equal(elw.ccurv, el.ccurv)


def test_writere2_3d():
    import pymech.neksuite as ns

    fin = "./tests/nek/box3d.re2"
    fout = "./test_2.re2"
    mesh = ns.readre2(fin)
    status = ns.writere2(fout, mesh)

    assert status == 0

    meshw = ns.readre2(fout)

    assert meshw.ndim == mesh.ndim
    assert meshw.nel == mesh.nel
    assert meshw.ncurv == mesh.ncurv
    assert meshw.nbc == mesh.nbc
    assert meshw.var == mesh.var
    assert meshw.lr1 == mesh.lr1
    assert meshw.wdsz == 8
    for (el, elw) in zip(mesh.elem, meshw.elem):
        npt.assert_array_equal(elw.pos, el.pos)
        npt.assert_array_equal(elw.bcs, el.bcs)
        npt.assert_array_equal(elw.curv, el.curv)
        npt.assert_array_equal(elw.ccurv, el.ccurv)


# ------------------------------------------------------------------------------
# test simson scripts
#
def test_readdns():
    import pymech.simsonsuite as ss

    fname = "./tests/simson/channel3D_t10000v.u"
    field = ss.readdns(fname)

    assert field.endian == "little"
    assert field.istep == []
    assert field.lr1 == [48, 65, 48]
    assert field.ndim == 3
    assert field.nel == 1
    assert field.var == [3, 3, 0, 0, 0]
    assert field.wdsz == 8
    assert (field.time - 10000.439742009798) < 1e-3


def test_readplane():
    import pymech.simsonsuite as ss

    fname = "./tests/simson/u.plane"
    x, d, nn, ndim = ss.readplane(fname)

    assert (x[0][1][0] - 0.06875) < 1e-3
    assert (d[0][1] - 0.0034688727137604305) < 1e-3
    assert nn[0] == 97.0
    assert nn[1] == 97.0
    assert ndim == 2


# ------------------------------------------------------------------------------
# test xarray dataset interface
#
def test_nekdataset():
    import pymech.dataset as pd

    fname = "./tests/nek/channel3D_0.f00001"
    ds = pd.open_dataset(fname)

    assert tuple(ds.dims.values()) == (64, 64, 64)
    assert math.isclose(ds.x.max(), 2 * math.pi, abs_tol=1e-6)
    assert math.isclose(ds.y.max(), 1.0, abs_tol=1e-6)
    assert math.isclose(ds.z.max(), math.pi, abs_tol=1e-6)
    assert math.isclose(ds.time, 0.2, abs_tol=1e-6)
