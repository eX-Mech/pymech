import math
import time
from textwrap import dedent

import numpy as np
from numpy import testing as npt
import pytest

from pymech.log import logger


# ------------------------------------------------------------------------------
# test nek scripts
#
def test_readnek(test_data_dir):
    import pymech.neksuite as ns

    fname = f"{test_data_dir}/nek/channel3D_0.f00001"
    field = ns.readnek(fname)

    assert field.endian == "little"
    assert field.istep == 10
    assert field.lr1 == (8, 8, 8)
    assert field.ndim == 3
    assert field.nel == 512
    assert field.var == (3, 3, 1, 0, 0)
    assert field.wdsz == 4
    assert (field.time - 0.2) < 1e-3
    representation = dedent(
        """\
    <pymech.core.HexaData>
    Dimensions:    3
    Precision:     4 bytes
    Mesh limits:
      * x:         (0.0, 6.2831854820251465)
      * y:         (-1.0, 1.0)
      * z:         (0.0, 3.1415927410125732)
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


def test_writenek(test_data_dir, tmpdir):
    import pymech.neksuite as ns

    fname = f"{test_data_dir}/nek/channel3D_0.f00001"
    time0 = time.perf_counter()
    field = ns.readnek(fname)
    time1 = time.perf_counter()

    fnamew = f"{tmpdir}/test_0.f00001"
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


def test_writenek_big_endian(test_data_dir, tmpdir):
    import pymech.neksuite as ns

    fname = f"{test_data_dir}/nek/channel3D_0.f00001"
    field = ns.readnek(fname)

    field.endian = "big"
    fnamew = f"{tmpdir}/test_0_big.f00001"
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


def test_readnek_scalars(test_data_dir):
    import pymech.neksuite as ns

    # 2D statistics file
    dtype = np.float32
    fname = f"{test_data_dir}/nek/stsabl0.f00001"
    field = ns.readnek(fname, dtype)

    ux_min, ux_max = field.lims.scal[0]
    assert math.isclose(ux_max, 5.3, abs_tol=0.1)
    assert field.elem[0].scal.dtype == np.dtype(dtype)


def test_writenek_scalars(test_data_dir, tmpdir):
    import pymech.neksuite as ns

    fname = f"{test_data_dir}/nek/stsabl0.f00001"
    field = ns.readnek(fname)

    fnamew = f"{tmpdir}/test_sts_0.f00001"
    status = ns.writenek(fnamew, field)

    assert status == 0

    fieldw = ns.readnek(fnamew)
    npt.assert_array_equal(field.lims.scal, fieldw.lims.scal)


def test_readnek_skip_vars(test_data_dir):
    import pymech.neksuite as ns

    fname = f"{test_data_dir}/nek/channel3D_0.f00001"

    field_all = ns.readnek(fname)
    field_skip_geom = ns.readnek(fname, skip_vars=("x", "y", "z"))
    field_skip_ux_uy = ns.readnek(fname, skip_vars=("ux", "uy"))

    for elem, elem_skip_geom, elem_skip_ux_uy in zip(
        field_all.elem, field_skip_geom.elem, field_skip_ux_uy.elem
    ):
        npt.assert_array_equal(elem.vel, elem_skip_geom.vel)
        npt.assert_array_equal(elem.pres, elem_skip_geom.pres)
        npt.assert_array_equal(elem.scal, elem_skip_geom.scal)
        npt.assert_array_equal(elem.pos, elem_skip_ux_uy.pos)
        npt.assert_array_equal(elem.pres, elem_skip_ux_uy.pres)
        npt.assert_array_equal(elem.scal, elem_skip_ux_uy.scal)
        with pytest.raises(AssertionError):
            npt.assert_array_equal(elem.pos, elem_skip_geom.pos)

        with pytest.raises(AssertionError):
            npt.assert_array_equal(elem.vel, elem_skip_ux_uy.vel)


def test_readnek_skip_scalars(test_data_dir):
    import pymech.neksuite as ns

    fname = f"{test_data_dir}/nek/stsabl0.f00001"
    field_all = ns.readnek(fname)
    field_skip = ns.readnek(fname, skip_vars=("s02, s04"))
    for elem, elem_skip in zip(field_all.elem, field_skip.elem):
        npt.assert_array_equal(elem.scal[0], elem_skip.scal[0])
        with pytest.raises(AssertionError):
            npt.assert_array_equal(elem.scal[1], elem_skip.scal[1])  # s02

        npt.assert_array_equal(elem.scal[2], elem_skip.scal[2])

        with pytest.raises(AssertionError):
            npt.assert_array_equal(elem.scal[3], elem_skip.scal[3])  # s04


def test_readrea(test_data_dir):
    import pymech.neksuite as ns

    fname = f"{test_data_dir}/nek/2D_section_R360.rea"
    field = ns.readrea(fname)

    assert field.lr1 == [2, 2, 1]
    assert field.ndim == 2
    assert field.nel == 1248
    assert abs(field.elem[0].pos[0][0][0][0] - 0.048383219999999998) < 1e-3
    assert abs(field.elem[887].curv[1, 0] - 1.21664) < 1e-3
    assert field.elem[887].ccurv[1] == "C"

    fname = f"{test_data_dir}/nek/m3j_bf_test.rea"
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


def test_writerea(test_data_dir, tmpdir):
    import pymech.neksuite as ns

    fname = f"{test_data_dir}/nek/2D_section_R360.rea"
    field = ns.readrea(fname)

    fnamew = f"{tmpdir}/test.rea"
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

    fname = f"{test_data_dir}/nek/m3j_bf_test.rea"
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


def test_merge(test_data_dir):
    import pymech.neksuite as ns
    import copy

    fname = f"{test_data_dir}/nek/box3d.rea"
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
        for iel, el in enumerate(mesh1.elem):
            for iface in range(6):
                assert el.bcs[ibc, iface][1] == iel + 1
                assert el.bcs[ibc, iface][2] == iface + 1


def test_readre2(test_data_dir):
    import pymech.neksuite as ns

    # The .re2 has been generated with reatore2 and contains the same data
    # except for the internal boundary conditions.
    # Assuming that `readrea` is correct, this checks id the .re2 file is read correctly too.
    #
    frea = f"{test_data_dir}/nek/2D_section_R360.rea"
    fre2 = f"{test_data_dir}/nek/2D_section_R360.re2"
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
    for el, elw in zip(meshrea.elem, meshre2.elem):
        npt.assert_allclose(elw.pos, el.pos)
        npt.assert_array_equal(elw.bcs, el.bcs)
        npt.assert_allclose(elw.curv, el.curv)
        npt.assert_array_equal(elw.ccurv, el.ccurv)


def test_readre2_3d(test_data_dir):
    import pymech.neksuite as ns

    # same test as test_readre2(), but with a 3D mesh.

    frea = f"{test_data_dir}/nek/box3d.rea"
    fre2 = f"{test_data_dir}/nek/box3d.re2"
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
    for el, elw in zip(meshrea.elem, meshre2.elem):
        npt.assert_allclose(elw.pos, el.pos)
        npt.assert_array_equal(elw.bcs, el.bcs)
        npt.assert_allclose(elw.curv, el.curv)
        npt.assert_array_equal(elw.ccurv, el.ccurv)


def test_writere2(test_data_dir, tmpdir):
    import pymech.neksuite as ns

    fin = f"{test_data_dir}/nek/2D_section_R360.re2"
    fout = f"{tmpdir}/test_1.re2"
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
    for el, elw in zip(mesh.elem, meshw.elem):
        npt.assert_array_equal(elw.pos, el.pos)
        npt.assert_array_equal(elw.bcs, el.bcs)
        npt.assert_array_equal(elw.curv, el.curv)
        npt.assert_array_equal(elw.ccurv, el.ccurv)


def test_writere2_3d(test_data_dir, tmpdir):
    import pymech.neksuite as ns

    fin = f"{test_data_dir}/nek/box3d.re2"
    fout = f"{tmpdir}/test_2.re2"
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
    for el, elw in zip(mesh.elem, meshw.elem):
        npt.assert_array_equal(elw.pos, el.pos)
        npt.assert_array_equal(elw.bcs, el.bcs)
        npt.assert_array_equal(elw.curv, el.curv)
        npt.assert_array_equal(elw.ccurv, el.ccurv)


def test_generate_internal_bcs(test_data_dir):
    import pymech.neksuite as ns
    import pymech.meshtools as mt

    # The rea and re2 meshes should be identical with the exception of internal boundary conditions.
    # The idea is to reconstruct the internal BCs of the re2 and compare with the .rea. They should be identical.
    frea = f"{test_data_dir}/nek/box3d.rea"
    fre2 = f"{test_data_dir}/nek/box3d.re2"
    meshrea = ns.readrea(frea)
    meshre2 = ns.readre2(fre2)
    nconnect = mt.generate_internal_bcs(meshre2)
    assert nconnect == 54  # This is a 3x3x3 box
    for ela, el2 in zip(meshrea.elem, meshre2.elem):
        npt.assert_array_equal(el2.bcs, ela.bcs)
    assert meshre2.check_connectivity()


def test_delete_internal_bcs(test_data_dir):
    import pymech.neksuite as ns
    import pymech.meshtools as mt

    # The rea and re2 meshes should be identical with the exception of internal boundary conditions.
    frea = f"{test_data_dir}/nek/box3d.rea"
    fre2 = f"{test_data_dir}/nek/box3d.re2"
    meshrea = ns.readrea(frea)
    meshre2 = ns.readre2(fre2)
    ndelete = mt.delete_internal_bcs(meshrea)
    assert (
        ndelete == 108
    )  # This is a 3x3x3 box, and each connection is deleted twice, one for each connected element
    for ela, el2 in zip(meshrea.elem, meshre2.elem):
        npt.assert_array_equal(el2.bcs, ela.bcs)
    assert meshrea.check_connectivity()


def test_extrude(test_data_dir):
    import pymech.neksuite as ns
    import pymech.meshtools as mt

    fname = f"{test_data_dir}/nek/2D_section_R360.re2"
    nz = 4
    z = np.linspace(-1, 1, nz + 1)
    mesh = ns.readre2(fname)
    mesh3D = mt.extrude(mesh, z)

    assert mesh3D.ndim == 3
    assert mesh3D.nel == mesh.nel * nz
    # curves duplicated on each side of each element
    assert mesh3D.ncurv == mesh.ncurv * nz * 2
    # check new periodic BCs in particular
    assert mesh3D.check_connectivity()


def test_extrude_refine(test_data_dir):
    import pymech.neksuite as ns
    import pymech.meshtools as mt
    from itertools import product

    fnameI = f"{test_data_dir}/nek/box2d.re2"
    mesh2D = ns.readre2(fnameI)
    mt.generate_internal_bcs(mesh2D)

    zmin = 0
    zmax = 6
    n = 16
    bc1 = ["P"]
    bc2 = ["P"]
    imesh_high = 0
    funpar = [0.5, 1.5]

    def fun_line(xpos, ypos, rlim):
        return ypos - rlim

    fun = [fun_line, fun_line]
    z = np.linspace(zmin, zmax, n + 1)

    # test both with and without internal connectivity
    mesh3D = mt.extrude_refine(
        mesh2D,
        z,
        bc1=bc1,
        bc2=bc2,
        fun=fun,
        funpar=funpar,
        imesh_high=imesh_high,
        internal_bcs=False,
    )
    assert mesh3D.check_connectivity()
    # check that we haven't introduced any dummy conditions
    for el, iface in product(mesh3D.elem, range(6)):
        assert el.bcs[0, iface][0] != "con"
        assert el.bcs[0, iface][0] != "E"

    mesh3D = mt.extrude_refine(
        mesh2D,
        z,
        bc1=bc1,
        bc2=bc2,
        fun=fun,
        funpar=funpar,
        imesh_high=imesh_high,
        internal_bcs=True,
    )
    for el, iface in product(mesh3D.elem, range(6)):
        assert el.bcs[0, iface][0] != "con"

    assert mesh3D.ndim == 3
    assert mesh3D.nel == 336
    assert mesh3D.check_connectivity()
    assert mesh3D.check_bcs_present()


def test_gen_circle(test_data_dir):
    import pymech.meshtools as mt

    # try a tiny mesh
    mesh = mt.gen_circle(1, 0.5, 2, 2)
    assert mesh.check_connectivity()
    assert mesh.check_bcs_present()
    assert mesh.nel == 20

    # one with ns > no
    mesh = mt.gen_circle(1, 0.5, 9, 2)
    assert mesh.check_connectivity()
    assert mesh.check_bcs_present()
    assert mesh.nel == 153

    # and a big one, without internal BCs
    mesh = mt.gen_circle(1, 0.1, 10, 200, internal_bcs=False)
    assert mesh.check_connectivity()
    assert mesh.nel == 8100


def test_cylinder(test_data_dir):
    import pymech.meshtools as mt
    import numpy as np

    # test making a cylinder using `gen_circle` and `extrude`
    # with a temperature field
    mesh = mt.gen_circle(1, 0.5, 9, 2, bc=["W", "I"])
    assert mesh.nbc == 2
    z = np.linspace(-1, 1, 5)
    # with default (periodic) boundary conditions
    mesh3D = mt.extrude(mesh, z)
    assert mesh3D.nbc == 2
    mesh3D.elem[0].bcs[1, 4][0] == "P"
    # with custom boundary conditions
    mesh3D = mt.extrude(mesh, z, bc1=["W", "t"], bc2=["W", "t"])
    assert mesh3D.elem[0].bcs[1, 4][0] == "t"


def test_gen_box(test_data_dir):
    import pymech.meshtools as mt

    # make a simple box and check that everything works fine
    box = mt.gen_box(4, 5, -1.0, 1.0, -1.0, 1.2)
    assert box.check_connectivity()
    assert box.nel == 20

    # try also with a temperature field and different boundary conditions
    box = mt.gen_box(
        4,
        5,
        -1.0,
        1.0,
        -1.0,
        1.0,
        var=[2, 2, 1, 1, 0],
        bcs_xmin=["v", "t"],
        bcs_xmax=["O", "I"],
        bcs_ymin=["W", "I"],
        bcs_ymax=["W", "I"],
    )
    assert box.check_connectivity()
    assert box.nbc == 2


def test_map2D(test_data_dir):
    import pymech.meshtools as mt

    # make a very twisted mesh and check that it ends up the way we want
    def mapping(x, y):
        return (
            2 * x + 0.2 * y + 0.1 * math.sin(math.pi * y),
            1.2 * y - 0.1 * x - 0.2 * math.sin(2 * math.pi * x),
        )

    n_centre = 5
    n_bl = 6
    s_param = 0.5
    radius = 1
    circle_mesh = mt.gen_circle(
        radius, s_param, n_centre, n_bl, var=[2, 2, 1, 1, 0], bc=["W", "I"]
    )
    circle_mesh = mt.map2D(circle_mesh, mapping)
    assert circle_mesh.check_connectivity()
    math.isclose(circle_mesh.elem[35].curv[0, 0], 1.4794589904304334)
    math.isclose(circle_mesh.elem[35].curv[3, 1], -0.2516125141385187)
    assert circle_mesh.elem[38].ccurv == [
        "m",
        "m",
        "m",
        "m",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
    ]

    # test with a box without any curvature
    # box of size 5×1
    xmin = 0
    xmax = 5
    ymin = 0
    ymax = 1
    # 20 × 10 elements resolution
    nx = 20
    ny = 10
    box = mt.gen_box(nx, ny, xmin, xmax, ymin, ymax)
    # change resolution such that the first element in y has a height of 0.025 instead of 0.1
    l0 = 0.025
    alpha = mt.exponential_refinement_parameter(l0, ymax, ny)

    def refinement_function(x, y):
        iy = ny * y / ymax
        return (x, l0 * iy**alpha)

    box = mt.map2D(box, refinement_function, curvature=False, boundary_curvature=False)
    box.check_connectivity()
    assert box.elem[0].ccurv == ["", "", "", "", "", "", "", "", "", "", "", ""]
    assert box.elem[25].ccurv == ["", "", "", "", "", "", "", "", "", "", "", ""]

    # and now with curvature activated on the boundary only
    box = mt.map2D(box, refinement_function, curvature=False)
    assert box.elem[0].ccurv == ["m", "", "", "m", "", "", "", "", "", "", "", ""]
    assert box.elem[25].ccurv == ["", "", "", "", "", "", "", "", "", "", "", ""]


# ------------------------------------------------------------------------------
# test simson scripts
#
def test_readdns(test_data_dir):
    import pymech.simsonsuite as ss

    fname = f"{test_data_dir}/simson/channel3D_t10000v.u"
    field = ss.readdns(fname)

    assert field.endian == "little"
    assert field.istep == []
    assert field.lr1 == [48, 65, 48]
    assert field.ndim == 3
    assert field.nel == 1
    assert field.var == [3, 3, 0, 0, 0]
    assert field.wdsz == 8
    assert (field.time - 10000.439742009798) < 1e-3


def test_readplane(test_data_dir):
    import pymech.simsonsuite as ss

    fname = f"{test_data_dir}/simson/u.plane"
    x, d, nn, ndim = ss.readplane(fname)

    assert (x[0][1][0] - 0.06875) < 1e-3
    assert (d[0][1] - 0.0034688727137604305) < 1e-3
    assert nn[0] == 97.0
    assert nn[1] == 97.0
    assert ndim == 2
