import copy
from itertools import product
from math import atan, cos, log, pi, sin, sqrt

import numpy as np

from pymech.core import HexaData
from pymech.log import logger


# ==============================================================================
def extrude(mesh: HexaData, z, bc1=None, bc2=None, internal_bcs=True):
    """Extrudes a 2D mesh into a 3D one

    Parameters
    ----------
    mesh : :class:`pymech.core.HexaData`
           2D mesh structure to extrude
    z :  float 1d array
           z coordinates at which to extrude the mesh, in increasing order
    bc1: str list
           A list of boundary conditions to use at the first end, one string per field. Defaults to periodic.
    bc2: str list
           A list of boundary conditions to use at the other end, one string per field. Defaults to periodic.
    internal_bcs : bool
           if True, build mesh connectivity using internal 'E' boundary conditions
           (note that those are not used by Nek5000 and will not be written to binary .re2 files).
    """

    # Set periodic boundary conditions by default if nothing else is requested
    if bc1 is None:
        bc1 = ["P"] * mesh.nbc
    if bc2 is None:
        bc2 = ["P"] * mesh.nbc

    if mesh.ndim != 2:
        raise ValueError("The mesh to extrude must be 2D")
    if mesh.lr1 != [2, 2, 1]:
        raise ValueError("Only mesh structures can be extruded (lr1 = [2, 2, 1])")
    if mesh.var[0] < 2:
        raise ValueError("The mesh to extrude must contain (x, y) geometry")
    # Is it possible to have periodic conditions for only some of the fields? If not, we should check for it too.
    for bc1_field, bc2_field in zip(bc1, bc2):
        if (bc1_field == "P" and bc2_field != "P") or (
            bc1_field != "P" and bc2_field == "P"
        ):
            raise ValueError(
                "Inconsistent boundary conditions: one end is periodic ('P') but the other isn't"
            )

    # copy the structure and make it 3D
    mesh3d = copy.deepcopy(mesh)
    if not internal_bcs:
        # if we don't build proper internal boundary conditions, we can get rid
        # of all of them to make the data cleaner
        delete_internal_bcs(mesh3d)
    mesh3d.lr1 = [2, 2, 2]
    mesh3d.var = [3, 0, 0, 0, 0]  # remove anything that isn't geometry for now
    nel2d = mesh.nel
    nz = len(z) - 1
    nel3d = mesh.nel * nz
    nbc = mesh.nbc
    mesh3d.nel = nel3d
    mesh3d.ndim = 3
    # The curved sides will also be extruded, one on each side of each element along nz
    mesh3d.ncurv = 2 * nz * mesh.ncurv

    # add extra copies of all elements
    for k in range(nz - 1):
        mesh_add = copy.deepcopy(mesh)
        # fix the indexing of periodic and internal conditions
        mesh_add.offset_connectivity((k + 1) * nel2d)
        mesh3d.elem = mesh3d.elem + mesh_add.elem

    # set the z locations and curvature
    for k in range(nz):
        for i in range(nel2d):
            iel = i + nel2d * k
            # replace the position arrays with 3D ones (with z empty)
            mesh3d.elem[iel].pos = np.zeros((3, 2, 2, 2))
            mesh3d.elem[iel].pos[:, :1, :, :] = mesh.elem[i].pos
            mesh3d.elem[iel].pos[:, 1:2, :, :] = mesh.elem[i].pos

            # fill in the z location
            z1 = z[k]
            z2 = z[k + 1]
            mesh3d.elem[iel].pos[2, 0, :, :] = z1
            mesh3d.elem[iel].pos[2, 1, :, :] = z2

            # extend curvature and correct it if necessary
            for icurv in range(4):
                curv_type = mesh3d.elem[iel].ccurv[icurv]
                # a 2D element has 4 edges that can be curved, numbered 0-3;
                # the extruded 3D element can have four more (on the other side), numbered 4-7
                mesh3d.elem[iel].ccurv[icurv + 4] = curv_type
                mesh3d.elem[iel].curv[icurv + 4] = mesh3d.elem[iel].curv[
                    icurv
                ]  # curvature params
                if curv_type == "m":
                    # in this case the midpoint is given. (x, y) is correct but z should be set to the proper value.
                    mesh3d.elem[iel].curv[icurv][2] = z1
                    mesh3d.elem[iel].curv[icurv + 4][2] = z2

    # fix the internal boundary conditions
    # the end boundary conditions will be overwritten later with the proper ones
    if internal_bcs:
        for iel, el in enumerate(mesh3d.elem):
            for ibc in range(nbc):
                el.bcs[ibc, 4][0] = "E"
                el.bcs[ibc, 4][1] = iel + 1
                el.bcs[ibc, 4][2] = 5
                el.bcs[ibc, 4][3] = iel - nel2d + 1
                el.bcs[ibc, 4][4] = 6
                el.bcs[ibc, 5][0] = "E"
                el.bcs[ibc, 5][1] = iel + 1
                el.bcs[ibc, 5][2] = 6
                el.bcs[ibc, 5][3] = iel + nel2d + 1
                el.bcs[ibc, 5][4] = 5
                # update the conditions for side faces
                for iface in range(4):
                    el.bcs[ibc, iface][1] = iel + 1
                    if el.bcs[ibc, iface][0] == "E":
                        # el.bcs[ibc, 0][1] ought to contain iel+1 once the mesh is valid
                        # but for now it should be off by a factor of nel2d because it is a copy of an element in the first slice
                        offset = iel - el.bcs[ibc, iface][1] + 1
                        el.bcs[ibc, iface][3] = el.bcs[ibc, iface][3] + offset

    # now fix the end boundary conditions
    # face 5 is at zmin and face 6 is at zmax (with Nek indexing, corresponding to 4 and 5 in Python)
    for i in range(nel2d):
        for ibc in range(nbc):
            i1 = i + (nz - 1) * nel2d  # index of the face on the zmax side
            mesh3d.elem[i].bcs[ibc, 4][0] = bc1[ibc]
            mesh3d.elem[i].bcs[ibc, 4][1] = i + 1
            mesh3d.elem[i].bcs[ibc, 4][2] = 5
            mesh3d.elem[i1].bcs[ibc, 5][0] = bc2[ibc]
            mesh3d.elem[i1].bcs[ibc, 5][1] = i1 + 1
            mesh3d.elem[i1].bcs[ibc, 5][2] = 6
            # fix the matching faces for the periodic conditions
            if bc1[ibc] == "P":
                mesh3d.elem[i].bcs[ibc, 4][3] = i1 + 1
                mesh3d.elem[i].bcs[ibc, 4][4] = 6
            else:
                mesh3d.elem[i].bcs[ibc, 4][3] = 0
                mesh3d.elem[i].bcs[ibc, 4][4] = 0
            if bc2[ibc] == "P":
                mesh3d.elem[i1].bcs[ibc, 5][3] = i + 1
                mesh3d.elem[i1].bcs[ibc, 5][4] = 5
            else:
                mesh3d.elem[i1].bcs[ibc, 5][3] = 0
                mesh3d.elem[i1].bcs[ibc, 5][4] = 0

    # return the extruded mesh
    return mesh3d


# ==============================================================================
def extrude_refine(
    mesh2D,
    z,
    bc1=None,
    bc2=None,
    fun=None,
    funpar=None,
    imesh_high=0,
    internal_bcs=True,
):
    r"""Extrudes a 2D mesh into a 3D one, following the pattern

    .. code-block::

         _____ _____ _____ _____
        |     |     |     |     |
        |     |     |     |     |
        |     |     |     |     |
        |_____|_____|_____|_____|
        |    /|\    |    /|\    |
        |__ / | \ __|__ / | \ __| (fun (with parameter funpar) should change change sign in the mid element)
        |  |  |  |  |  |  |  |  | (half of the mid elements are also divided in 2 in (x,y)-plane)
        |__|__|__|__|__|__|__|__|
        |  |  |  |  |  |  |  |  |
        |  |  |  |  |  |  |  |  |
        |  |  |  |  |  |  |  |  | (imesh_high is the index of mesh with higher intended discretization in z)
        |__|__|__|__|__|__|__|__|

    The pattern is similar to "Picture Frame" of an ancient NEKTON manual (https://www.mcs.anl.gov/~fischer/Nek5000/nekmanual.pdf).
    If the mid elements have curvature, the extrusion might modify it. Do not split in regions where the value of curvature parameters is very important.

    Parameters
    ----------
    mesh2D : :class:`pymech.core.HexaData`
           2D mesh structure to extrude
    z : float array
        list of z values of the  most refined zones of the extruded mesh
    bc1: str list
           A list of boundary conditions to use at the first end, one string per field. Defaults to periodic.
    bc2: str list
           A list of boundary conditions to use at the other end, one string per field. Defaults to periodic.
    fun: function
         list of functions that define the splitting lines for different discretization meshes (default: empty, in which case the simple extrusion function `extrude` is called instead)
    funpar: list
          list of parameters for functions that define the splitting lines for different discretization meshes (default: None, equivalent to an array of zeroes)
    imesh_high : int
                 index of fun that defines the mesh with higher discretization. Example: 0, is the most internal mesh; 1 is the second most internal mesh, etc (default: the most internal mesh, imesh_high=0)
    """

    # Set periodic boundary conditions by default if nothing else is requested
    if bc1 is None:
        bc1 = ["P"] * mesh2D.nbc
    if bc2 is None:
        bc2 = ["P"] * mesh2D.nbc

    # Consistency checks: Initial grid
    if mesh2D.ndim != 2:
        raise ValueError("The mesh to extrude must be 2D")
    if mesh2D.lr1 != [2, 2, 1]:
        raise ValueError("Only mesh structures can be extruded (lr1 = [2, 2, 1])")
    if mesh2D.var[0] < 2:
        raise ValueError("The mesh to extrude must contain (x, y) geometry")
    # Consistency checks: Periodic boundary condition
    for bc1_field, bc2_field in zip(bc1, bc2):
        if (bc1_field == "P" and bc2_field != "P") or (
            bc1_field != "P" and bc2_field == "P"
        ):
            raise ValueError(
                "Inconsistent boundary conditions: one end is 'P' but the other isn't"
            )

    # Consistency checks: Functions that define the splitting lines
    nsplit = len(fun)
    if funpar is not None and len(funpar) != nsplit:
        raise ValueError(
            f"The length of funpar ({len(funpar)}) must match the length of par ({nsplit})!"
        )

    # number of elements in the z direction
    nz = len(z) - 1

    # Consistency checks: if nz is divided by 4 (or 8, or 16, etc)
    if (nz % 2 ** abs(imesh_high + 1) != 0) or (
        nz % 2 ** abs(nsplit - imesh_high + 1) != 0
    ):
        raise ValueError(
            f"Inconsistent elements to extrude: the number of elements ({nz}) must be a multiple of {max([2**abs(imesh_high + 1), 2**abs(nsplit - imesh_high + 1)])}"
        )

    # If fun is not defined, there is no splitting to be done. Call simple extrusion and end routine
    if fun is None:
        logger.info("Splitting function not defined. Calling simple extrusion routine.")
        mesh3D = extrude(mesh2D, z, bc1, bc2)
        return mesh3D

    mesh2D_ext = copy.deepcopy(mesh2D)
    # get rid of internal boundary conditions if they are not requested
    if not internal_bcs:
        delete_internal_bcs(mesh2D_ext)

    meshes2D = []  # list of 2D meshes
    meshes3D = []  # list of 3D meshes

    # Splitting 2D meshes
    # sorts element of the 2D mesh into layers: most refined, first split, second most refined,
    # second split, etc, until the least refined stored into mesh2D_ext.
    for k in range(nsplit):
        meshes2D.append(copy.deepcopy(mesh2D_ext))
        meshes2D.append(copy.deepcopy(mesh2D_ext))

        elems_int = []
        elems_mid = []
        elems_ext = []

        for iel in range(mesh2D_ext.nel):
            it = 0
            xvec = np.zeros((4, 1))
            yvec = np.zeros((4, 1))
            rvec = np.zeros((4, 1))
            for i in range(2):
                for j in range(2):
                    xvec[it] = mesh2D_ext.elem[iel].pos[0, 0, j, i]
                    yvec[it] = mesh2D_ext.elem[iel].pos[1, 0, j, i]
                    rvec[it] = fun[k](xvec[it], yvec[it], funpar[k])
                    it += 1
            if (
                max(rvec) <= 0.0
            ):  # the element belongs to the internal (more refined) mesh
                elems_int.append(iel)
            elif (
                min(rvec) > 0.0
            ):  # the element belongs to the external (less refined) mesh
                elems_ext.append(iel)
            else:  # the element belongs to the intermediate mesh and will be split
                elems_mid.append(iel)

        external_bc = ""
        if internal_bcs:
            # dummy "external" boundary condition that will be used between parts of meshes
            # to be merged together later in order to mark the faces to merge.
            # This is only necessary if we want to build the internal connectivity.
            external_bc = "con"
        keep_elements(meshes2D[2 * k], elems_int, external_bc=external_bc)
        keep_elements(meshes2D[2 * k + 1], elems_mid, external_bc=external_bc)
        keep_elements(mesh2D_ext, elems_ext, external_bc=external_bc)
        logger.debug(f"Mesh2D {2 * k} elements: {meshes2D[2 * k].nel}")
        logger.debug(f"Mesh2D {2 * k + 1} elements: {meshes2D[2 * k + 1].nel}")

    # End of splitting, remaining is the last mesh: Mesh_ext
    logger.debug(f"Mesh2Dext elements: {mesh2D_ext.nel}")

    # Update curvature metadata for all sub-meshes
    mesh2D_ext.update_ncurv()
    for mesh_part in meshes2D:
        mesh_part.update_ncurv()

    # Extruding meshes
    logger.info("Extruding meshes")
    for k in range(nsplit):
        # divide number of elements by 2**(coarsening level)
        n_local = nz // 2 ** abs(k - imesh_high)
        # select z coordinates for coarsened elements
        z_local = z[:: int(2 ** abs(k - imesh_high))]

        if k < imesh_high:

            def fun_local(xpos, ypos, rlim):
                return -fun[k](xpos, ypos, rlim)

            n_mid = 2 * n_local
            z_mid = z[:: int(2 ** abs(k - imesh_high + 1))]
        else:
            fun_local = fun[k]
            n_mid = n_local
            z_mid = z_local

        if n_mid % 4 != 0:
            raise ValueError(
                f"Inconsistent elements to extrude: n ({n_mid}) is not a multiple of 4."
            )

        meshes3D.append(
            extrude(
                meshes2D[2 * k], z_local, bc1=bc1, bc2=bc2, internal_bcs=internal_bcs
            )
        )
        meshes3D.append(
            extrude_mid(
                meshes2D[2 * k + 1],
                z_mid,
                bc1,
                bc2,
                fun_local,
                funpar=funpar[k],
                internal_bcs=internal_bcs,
            )
        )

        logger.debug(f"Mesh3D {2 * k} elements: {meshes3D[2 * k].nel}")
        logger.debug(f"Mesh3D {2 * k + 1} elements: {meshes3D[2 * k + 1].nel}")

    n_local = nz // 2 ** abs(nsplit - imesh_high)
    z_local = z[:: int(2 ** abs(nsplit - imesh_high))]
    mesh3D_ext = extrude(
        mesh2D_ext, z_local, bc1=bc1, bc2=bc2, internal_bcs=internal_bcs
    )

    logger.debug(f"Mesh3Dext elements: {mesh3D_ext.nel}")

    # Merging meshes
    logger.info("Merging meshes")
    mesh3D = mesh3D_ext
    for mesh_part in meshes3D:
        mesh3D.merge(mesh_part, ignore_all_bcs=not internal_bcs)
    logger.info(f"Merging done. Total elements: {mesh3D.nel}")

    # update curve metadata
    mesh3D.update_ncurv()

    # return the extruded mesh
    return mesh3D


# =================================================================================


def extrude_mid(mesh, z, bc1, bc2, fun, funpar=0.0, internal_bcs=True):
    r"""Extrudes the mid elments of the 2D mesh into a 3D one. Following the pattern:

    .. code-block::

         _____ _____ _____ _____
        |1   /|\   4|    /|\    |
        |__ / | \ __|__ / | \ __| (fun (with parameter funpar) should change change sign in the mid element)
        |0 |2 |3 | 5|  |  |  |  | (half of the mid elements are also divided in 2 in (x, y)-plane)
        |__|__|__|__|__|__|__|__| (numbers in the figure indicate the indices (iel + 0; iel + 1; etc))

    Parameters
    ----------
    mesh : :class:`pymech.core.HexaData`
           2D mesh structure to extrude
    z : float
        list of z values of the nodes of the elements of the extruded mesh in the high discretization region (len(z)-1 must be divide by 4)
    bc1: str list
           A list of boundary conditions to use at the first end, one string per field.
    bc2: str list
           A list of boundary conditions to use at the other end, one string per field.
    fun : function
          function that define the splitting lines for different discretization meshes
    funpar : not defined, depends on the function
             parameter for functions that define the splitting lines for different discretization meshes (default: zero, can be used for when funpar is not needed inside fun)

    Suggestion: see function extrude_split to understand how to call extrude_mid
    """

    # Consistency checks: Initial grid
    if mesh.ndim != 2:
        raise ValueError("The mesh to extrude must be 2D")
    if mesh.lr1 != [2, 2, 1]:
        raise ValueError("Only mesh structures can be extruded (lr1 = [2, 2, 1])")
    if mesh.var[0] < 2:
        raise ValueError("The mesh to extrude must contain (x, y) geometry")
    for bc1_field, bc2_field in zip(bc1, bc2):
        if (bc1_field == "P" and bc2_field != "P") or (
            bc1_field != "P" and bc2_field == "P"
        ):
            raise ValueError(
                "Inconsistent boundary conditions: one end is periodic ('P') but the other isn't"
            )

    nz = len(z) - 1
    z1 = np.zeros((nz, 1))
    z2 = np.zeros((nz, 1))
    z1 = z[0:nz]
    z2 = z[1 : nz + 1]

    if nz % 4 != 0:
        raise ValueError("Inconsistent elements to extrude: nz must be divided by 4")

    # copy the structure and make it 3D
    mesh3d = copy.deepcopy(mesh)

    # add extra copies of all elements
    for k in range(6 * (nz // 4) - 1):
        mesh_add = copy.deepcopy(mesh)
        mesh3d.merge(mesh_add, ignore_all_bcs=True)

    # make the mesh 3D
    mesh3d.lr1 = [2, 2, 2]
    mesh3d.var = [3, 0, 0, 0, 0]  # remove anything that isn't geometry for now
    nel2d = mesh.nel
    nel3d = (
        nel2d * 6 * (nz // 4)
    )  # every mid-extrusion of a 2d-element creates 6 3d-elements, while the usual extrusion creates 4 in the high discretized grid and 2 in the low discretized grid
    nbc = mesh.nbc
    mesh3d.nel = nel3d
    mesh3d.ndim = 3
    # The curved sides will also be extruded, one on each side of each element along nz
    mesh3d.ncurv = 2 * nz * mesh.ncurv

    # set the x, y, z locations and curvature
    for k in range(0, nz, 4):
        for i in range(nel2d):
            iel = 6 * (i + nel2d * (k // 4))
            for iell in range(6):
                mesh3d.elem[iel + iell] = copy.deepcopy(mesh.elem[i])
                mesh3d.elem[iel + iell].pos = np.zeros((3, 2, 2, 2))

            xvec = np.zeros((2, 2))
            yvec = np.zeros((2, 2))
            rvec = np.zeros((2, 2))
            index_lo = np.zeros(
                (2, 2), dtype=int
            )  # index [ix, iy] of the two low points
            index_hi = np.zeros(
                (2, 2), dtype=int
            )  # index [ix, iy] of the two high points
            iindex_lo = 0  # index of low points (among the four corners), that have a negative image by the function
            iindex_hi = (
                0  # index of high points, that have a positive image by the function
            )
            iedgelat = np.zeros((2), dtype=int)
            iedgeconlat = np.zeros((2), dtype=int)
            # find which of the four corners are low points and high points; there must be two of each.
            for ii in range(2):
                for jj in range(2):
                    xvec[jj, ii] = mesh.elem[i].pos[0, 0, jj, ii]
                    yvec[jj, ii] = mesh.elem[i].pos[1, 0, jj, ii]
                    rvec[jj, ii] = fun(xvec[jj, ii], yvec[jj, ii], funpar)
                    if rvec[jj, ii] <= 0.0:
                        if iindex_lo > 1:
                            raise ValueError(
                                "Mid element not consistent. Criteria must divide elements with 2 points on each side."
                            )
                        index_lo[iindex_lo, :] = [jj, ii]
                        iindex_lo += 1
                    else:
                        if iindex_hi > 1:
                            raise ValueError(
                                "Mid element not consistent. Criteria must divide elements with 2 points on each side."
                            )
                        index_hi[iindex_hi, :] = [jj, ii]
                        iindex_hi += 1
            if (iindex_lo != 2) or (iindex_hi != 2):
                raise ValueError(
                    "Mid element not consistent. Criteria must divide elements with 2 points on each side."
                )

            # find the indices of edges, for curvature and boundary condition
            #
            # iedgehi is the index of the edge of element iel + 0 that is the intersection between iel + 0 and iel + 1 (high edge).
            # iedgelo is the index of the edge of element iel + 1 that is the intersection (low edge).
            # iedgelat are the indices of the lateral (splitted) edges.
            # iedgeconlat are the indices of the edges (edge in z-direction) of elements iel + 2
            # and iel + 3 that connect to the respective lateral edges.
            if (index_lo[0, :] == [0, 0]).all():
                if (index_hi[0, :] == [0, 1]).all():
                    iedgehi = 1
                    iedgelo = 3
                    iedgelat[0] = 0
                    iedgelat[1] = 2
                    iedgeconlat[0] = 9
                    iedgeconlat[1] = 10
                else:
                    iedgehi = 2
                    iedgelo = 0
                    iedgelat[0] = 3
                    iedgelat[1] = 1
                    iedgeconlat[0] = 11
                    iedgeconlat[1] = 10
            elif (index_lo[0, :] == [1, 0]).all():
                iedgehi = 0
                iedgelo = 2
                iedgelat[0] = 3
                iedgelat[1] = 1
                iedgeconlat[0] = 8
                iedgeconlat[1] = 9
            elif (index_lo[0, :] == [0, 1]).all():
                iedgehi = 3
                iedgelo = 1
                iedgelat[0] = 0
                iedgelat[1] = 2
                iedgeconlat[0] = 8
                iedgeconlat[1] = 11

            # find x and y locations
            poslo = mesh.elem[i].pos[:, :1, index_lo[:, 0], index_lo[:, 1]]
            poshi = mesh.elem[i].pos[:, :1, index_hi[:, 0], index_hi[:, 1]]
            # mid position is influenced by curvature
            posmid = 0.5 * (
                mesh.elem[i].pos[:, :1, index_lo[:, 0], index_lo[:, 1]]
                + mesh.elem[i].pos[:, :1, index_hi[:, 0], index_hi[:, 1]]
            )
            for ilat in range(2):
                # finds the mid points of lateral edges (also considering curvature)
                posmid[:, 0, ilat] = edge_mid(mesh.elem[i], iedgelat[ilat])

            # fill in the x and y location
            mesh3d.elem[iel].pos[:, 0:2, index_lo[:, 0], index_lo[:, 1]] = poslo
            mesh3d.elem[iel].pos[:, 0:2, index_hi[:, 0], index_hi[:, 1]] = posmid

            mesh3d.elem[iel + 1].pos[:, 0:2, index_lo[:, 0], index_lo[:, 1]] = posmid
            mesh3d.elem[iel + 1].pos[:, 0:2, index_hi[:, 0], index_hi[:, 1]] = poshi

            mesh3d.elem[iel + 2].pos[:, 0:2, index_lo[:, 0], index_lo[:, 1]] = poslo
            mesh3d.elem[iel + 2].pos[:, :1, index_hi[:, 0], index_hi[:, 1]] = posmid
            mesh3d.elem[iel + 2].pos[:, 1:2, index_hi[:, 0], index_hi[:, 1]] = poshi

            mesh3d.elem[iel + 3].pos[:, 0:2, index_lo[:, 0], index_lo[:, 1]] = poslo
            mesh3d.elem[iel + 3].pos[:, :1, index_hi[:, 0], index_hi[:, 1]] = poshi
            mesh3d.elem[iel + 3].pos[:, 1:2, index_hi[:, 0], index_hi[:, 1]] = posmid

            mesh3d.elem[iel + 4].pos[:, 0:2, index_lo[:, 0], index_lo[:, 1]] = posmid
            mesh3d.elem[iel + 4].pos[:, 0:2, index_hi[:, 0], index_hi[:, 1]] = poshi

            mesh3d.elem[iel + 5].pos[:, 0:2, index_lo[:, 0], index_lo[:, 1]] = poslo
            mesh3d.elem[iel + 5].pos[:, 0:2, index_hi[:, 0], index_hi[:, 1]] = posmid

            # fill in the z location
            mesh3d.elem[iel].pos[2, 0, :, :] = z1[k]
            mesh3d.elem[iel].pos[2, 1, :, :] = z2[k]

            mesh3d.elem[iel + 1].pos[2, 0, :, :] = z1[k]
            mesh3d.elem[iel + 1].pos[2, 1, index_lo[0, 0], index_lo[0, 1]] = z2[k]
            mesh3d.elem[iel + 1].pos[2, 1, index_lo[1, 0], index_lo[1, 1]] = z2[k]
            mesh3d.elem[iel + 1].pos[2, 1, index_hi[0, 0], index_hi[0, 1]] = z2[k + 1]
            mesh3d.elem[iel + 1].pos[2, 1, index_hi[1, 0], index_hi[1, 1]] = z2[k + 1]

            mesh3d.elem[iel + 2].pos[2, 0, :, :] = z1[k + 1]
            mesh3d.elem[iel + 2].pos[2, 1, :, :] = z2[k + 1]

            mesh3d.elem[iel + 3].pos[2, 0, :, :] = z1[k + 2]
            mesh3d.elem[iel + 3].pos[2, 1, :, :] = z2[k + 2]

            mesh3d.elem[iel + 4].pos[2, 1, :, :] = z2[k + 3]
            mesh3d.elem[iel + 4].pos[2, 0, index_lo[0, 0], index_lo[0, 1]] = z1[k + 3]
            mesh3d.elem[iel + 4].pos[2, 0, index_lo[1, 0], index_lo[1, 1]] = z1[k + 3]
            mesh3d.elem[iel + 4].pos[2, 0, index_hi[0, 0], index_hi[0, 1]] = z1[k + 2]
            mesh3d.elem[iel + 4].pos[2, 0, index_hi[1, 0], index_hi[1, 1]] = z1[k + 2]

            mesh3d.elem[iel + 5].pos[2, 0, :, :] = z1[k + 3]
            mesh3d.elem[iel + 5].pos[2, 1, :, :] = z2[k + 3]

            for icurv in range(4):
                # Curvature type 's' acts on faces, not edges. Does not make sense for extruded mesh. Changing to 'C'.
                if mesh.elem[i].ccurv[icurv] == "s":
                    logger.warn(
                        f"Curvature s on element {i + 1}. Not consistent with extrusion, changing to C"
                    )
                    mesh.elem[i].ccurv[icurv] == "C"
                    mesh.elem[i].ccurv[icurv][0] = mesh.elem[i].ccurv[icurv][4]
                    mesh.elem[i].ccurv[icurv][1:4] = 0.0
                elif (
                    (mesh.elem[i].ccurv[icurv] != "")
                    and (mesh.elem[i].ccurv[icurv] != "m")
                    and (mesh.elem[i].ccurv[icurv] != "C")
                ):
                    logger.warn(f"Curvature unknown on element {i + 1}")

            # extend curvature and correct it if necessary
            # calculate coordinates of midsize-node for every curvature type, except both empty (even if both are 'C'). Then find radius, if applicable. 'm' takes precendence over 'C'.
            # if mesh.elem[i].ccurv[iedgehi] != mesh.elem[i].ccurv[iedgelo]:
            if mesh.elem[i].ccurv[iedgehi] != "" or mesh.elem[i].ccurv[iedgelo] != "":
                midpointhi = edge_mid(mesh.elem[i], iedgehi)
                midpointlo = edge_mid(mesh.elem[i], iedgelo)
                midpointmid = 0.5 * (posmid[:, 0, 0] + posmid[:, 0, 1]) + 0.5 * (
                    midpointhi
                    - 0.5 * (poshi[:, 0, 0] + poshi[:, 0, 1])
                    + midpointlo
                    - 0.5 * (poslo[:, 0, 0] + poslo[:, 0, 1])
                )
                if (
                    mesh.elem[i].ccurv[iedgehi] == "m"
                    or mesh.elem[i].ccurv[iedgelo] == "m"
                ):
                    mesh3d.elem[iel].ccurv[iedgehi] = "m"
                    mesh3d.elem[iel + 1].ccurv[iedgelo] = "m"
                    mesh3d.elem[iel].curv[iedgehi][:3] = midpointmid
                    mesh3d.elem[iel + 1].curv[iedgelo][:3] = midpointmid
                elif (
                    mesh.elem[i].ccurv[iedgehi] == "C"
                    or mesh.elem[i].ccurv[iedgelo] == "C"
                ):
                    midpointmid[2] = z1[k]
                    curvmid = edge_circle(mesh3d.elem[iel], iedgehi, midpointmid)
                    if curvmid[0] == 0.0:
                        mesh3d.elem[iel].ccurv[iedgehi] = ""
                        mesh3d.elem[iel + 1].ccurv[iedgelo] = ""
                        mesh3d.elem[iel].curv[iedgehi][:4] = 0.0
                        mesh3d.elem[iel + 1].curv[iedgelo][:4] = 0.0
                    else:
                        curvmid[1:4] = 0.0
                        mesh3d.elem[iel].ccurv[iedgehi] = "C"
                        mesh3d.elem[iel + 1].ccurv[iedgelo] = "C"
                        mesh3d.elem[iel].curv[iedgehi][:4] = curvmid
                        mesh3d.elem[iel + 1].curv[iedgelo][:4] = -curvmid
                else:
                    logger.warn(
                        f"Splitted element curvature unknown on element {i + 1} of mid mesh. Removing curvature."
                    )
                    # For cases not implemented, remove curvature
                    mesh3d.elem[iel].ccurv[iedgehi] = ""
                    mesh3d.elem[iel + 1].ccurv[iedgelo] = ""
                    mesh3d.elem[iel].curv[iedgehi] = 0.0 * mesh.elem[i].curv[iedgehi]
                    mesh3d.elem[iel + 1].curv[iedgelo] = mesh3d.elem[iel].curv[iedgehi]

            for ilat in range(2):
                # Fixing curvature of edges divided in half. For curv_type == 'C', it is not a true extrusion - 'diagonal edges' not consistent with 'C'.
                if mesh.elem[i].ccurv[iedgelat[ilat]] == "m":
                    # coordinates of midsize-node is approximated. Considering an edge aligned with the x-axis, the position would be (ym_new = 3/4*ym_old = 1/2*ym_old(mean y-position) + 1/4*ym_old(distance)) at xm_new=(x2-x1)/4. This comes from parabolic curve (however, it is not exactly midsize)
                    dmid = (
                        (poshi[0, 0, ilat] - poslo[0, 0, ilat])
                        * (poslo[1, 0, ilat] - posmid[1, 0, ilat])
                        - (poslo[0, 0, ilat] - posmid[0, 0, ilat])
                        * (poshi[1, 0, ilat] - poslo[1, 0, ilat])
                    ) / (
                        (poshi[0, 0, ilat] - poslo[0, 0, ilat]) ** 2
                        + (poshi[1, 0, ilat] - poslo[1, 0, ilat]) ** 2
                    )
                    mesh3d.elem[iel].curv[iedgelat[ilat]][0] = 0.5 * (
                        poslo[0, 0, ilat] + posmid[0, 0, ilat]
                    ) + dmid / 4.0 * (poshi[1, 0, ilat] - poslo[1, 0, ilat])
                    mesh3d.elem[iel + 1].curv[iedgelat[ilat]][0] = 0.5 * (
                        posmid[0, 0, ilat] + poshi[0, 0, ilat]
                    ) + dmid / 4.0 * (poshi[1, 0, ilat] - poslo[1, 0, ilat])
                    mesh3d.elem[iel].curv[iedgelat[ilat]][1] = 0.5 * (
                        poslo[1, 0, ilat] + posmid[1, 0, ilat]
                    ) - dmid / 4.0 * (poshi[0, 0, ilat] - poslo[0, 0, ilat])
                    mesh3d.elem[iel + 1].curv[iedgelat[ilat]][1] = 0.5 * (
                        posmid[1, 0, ilat] + poshi[1, 0, ilat]
                    ) - dmid / 4.0 * (poshi[0, 0, ilat] - poslo[0, 0, ilat])
                elif mesh.elem[i].ccurv[iedgelat[ilat]] == "C":
                    # if the lateral edge has curvature 'C', the diagonal edges connected to it (iedgeconlat of elements iel + 2 and iel + 3) would have curvature 'C', which are inconsistent with edges 8-12 inside Nek5000 (because these are supposed to be edges in z-direction). Changing to curvature 'm' for elements iel + 1 (and iel + 2, iel + 3, iel + 2, iel + 4)
                    midpointlathi = edge_mid(mesh3d.elem[iel + 1], iedgelat[ilat])
                    mesh3d.elem[iel + 1].curv[iedgelat[ilat]][:3] = midpointlathi
                    mesh3d.elem[iel + 1].ccurv[iedgelat[ilat]] = "m"

            for icurv in range(4):
                # a 2D element has 4 edges that can be curved, numbered 0-3;
                # the extruded 3D element can have four more (on the other side), numbered 4-7
                for iell in range(6):
                    mesh3d.elem[iel + iell].ccurv[icurv + 4] = mesh3d.elem[
                        iel + iell
                    ].ccurv[icurv]
                    mesh3d.elem[iel + iell].curv[icurv + 4] = mesh3d.elem[
                        iel + iell
                    ].curv[
                        icurv
                    ]  # curvature params

            mesh3d.elem[iel + 2].curv[0:4] = mesh3d.elem[iel].curv[0:4]
            mesh3d.elem[iel + 3].curv[4:8] = mesh3d.elem[iel].curv[0:4]
            mesh3d.elem[iel + 5].curv[:] = mesh3d.elem[iel].curv
            mesh3d.elem[iel + 4].curv[:] = mesh3d.elem[iel + 1].curv
            mesh3d.elem[iel + 2].ccurv[0:4] = mesh3d.elem[iel].ccurv[0:4]
            mesh3d.elem[iel + 3].ccurv[4:8] = mesh3d.elem[iel].ccurv[0:4]
            mesh3d.elem[iel + 5].ccurv[:] = mesh3d.elem[iel].ccurv
            mesh3d.elem[iel + 4].ccurv[:] = mesh3d.elem[iel + 1].ccurv

            for icurv in range(4):
                # z should be set to the proper value.
                curv_type = mesh3d.elem[iel].ccurv[icurv]
                if curv_type == "m":
                    izcurv = 2
                    mesh3d.elem[iel].curv[icurv][izcurv] = z1[k]
                    mesh3d.elem[iel].curv[icurv + 4][izcurv] = z2[k]
                    mesh3d.elem[iel + 2].curv[icurv][izcurv] = z1[k + 1]
                    mesh3d.elem[iel + 2].curv[icurv + 4][izcurv] = z2[k + 1]
                    mesh3d.elem[iel + 3].curv[icurv][izcurv] = z1[k + 2]
                    mesh3d.elem[iel + 3].curv[icurv + 4][izcurv] = z2[k + 2]
                    mesh3d.elem[iel + 5].curv[icurv][izcurv] = z1[k + 3]
                    mesh3d.elem[iel + 5].curv[icurv + 4][izcurv] = z2[k + 3]

                curv_type = mesh3d.elem[iel + 1].ccurv[icurv]
                # curvature of iel + 1 may be different from iel because of diagonal edges
                if curv_type == "m":
                    izcurv = 2
                    mesh3d.elem[iel + 1].curv[icurv][izcurv] = z1[k]
                    mesh3d.elem[iel + 4].curv[icurv + 4][izcurv] = z2[k + 3]
                    if icurv == iedgehi:
                        mesh3d.elem[iel + 1].curv[iedgehi + 4][izcurv] = z2[k + 1]
                        mesh3d.elem[iel + 4].curv[iedgehi][izcurv] = z2[k + 1]
                    elif icurv == iedgelo:
                        mesh3d.elem[iel + 1].curv[iedgelo + 4][izcurv] = z1[k + 1]
                        mesh3d.elem[iel + 4].curv[iedgelo][izcurv] = z2[k + 2]
                    elif icurv == iedgelat[0] or icurv == iedgelat[1]:
                        mesh3d.elem[iel + 1].curv[icurv + 4][izcurv] = 0.5 * (
                            z1[k + 1] + z2[k + 1]
                        )
                        mesh3d.elem[iel + 4].curv[icurv][izcurv] = 0.5 * (
                            z2[k + 1] + z2[k + 2]
                        )

            # Fixing the curvature of 3d-edges in z-direction that connects to lateral edges in trapezoidal elements (all other edges in z-direction - indices 8 to 11 - do not have curvature)
            for ilat in range(2):
                mesh3d.elem[iel + 2].curv[iedgeconlat[ilat]] = mesh3d.elem[
                    iel + 1
                ].curv[iedgelat[ilat] + 4]
                mesh3d.elem[iel + 2].ccurv[iedgeconlat[ilat]] = mesh3d.elem[
                    iel + 1
                ].ccurv[iedgelat[ilat] + 4]
                mesh3d.elem[iel + 3].curv[iedgeconlat[ilat]] = mesh3d.elem[
                    iel + 4
                ].curv[iedgelat[ilat]]
                mesh3d.elem[iel + 3].ccurv[iedgeconlat[ilat]] = mesh3d.elem[
                    iel + 4
                ].ccurv[iedgelat[ilat]]

            # fix the internal boundary conditions
            # the end boundary conditions will be overwritten later with the proper ones
            for ibc in range(nbc):
                # set the conditions in faces normal to z
                if internal_bcs:
                    for iell in range(6):
                        mesh3d.elem[iel + iell].bcs[ibc, 4][0] = "E"
                        mesh3d.elem[iel + iell].bcs[ibc, 4][1] = iel + iell + 1
                        mesh3d.elem[iel + iell].bcs[ibc, 4][2] = 5
                        mesh3d.elem[iel + iell].bcs[ibc, 4][4] = 6
                        mesh3d.elem[iel + iell].bcs[ibc, 5][0] = "E"
                        mesh3d.elem[iel + iell].bcs[ibc, 5][1] = iel + iell + 1
                        mesh3d.elem[iel + iell].bcs[ibc, 5][2] = 6
                        mesh3d.elem[iel + iell].bcs[ibc, 5][4] = 5

                    mesh3d.elem[iel].bcs[ibc, 4][3] = iel + 1 + 5 - 6 * nel2d
                    mesh3d.elem[iel].bcs[ibc, 5][3] = iel + 1 + 2
                    mesh3d.elem[iel + 1].bcs[ibc, 4][3] = iel + 1 + 4 - 6 * nel2d
                    mesh3d.elem[iel + 2].bcs[ibc, 4][3] = iel + 1
                    mesh3d.elem[iel + 2].bcs[ibc, 5][3] = iel + 1 + 3
                    mesh3d.elem[iel + 3].bcs[ibc, 4][3] = iel + 1 + 2
                    mesh3d.elem[iel + 3].bcs[ibc, 5][3] = iel + 1 + 5
                    mesh3d.elem[iel + 4].bcs[ibc, 5][3] = iel + 1 + 1 + 6 * nel2d
                    mesh3d.elem[iel + 5].bcs[ibc, 4][3] = iel + 1 + 3
                    mesh3d.elem[iel + 5].bcs[ibc, 5][3] = iel + 1 + 6 * nel2d

                    # Correct internal bc for mid faces of elements.
                    mesh3d.elem[iel].bcs[ibc, iedgehi][0] = "E"
                    mesh3d.elem[iel].bcs[ibc, iedgehi][3] = iel + 1 + 1
                    mesh3d.elem[iel].bcs[ibc, iedgehi][4] = iedgelo + 1
                    mesh3d.elem[iel + 1].bcs[ibc, iedgelo][0] = "E"
                    mesh3d.elem[iel + 1].bcs[ibc, iedgelo][3] = iel + 1
                    mesh3d.elem[iel + 1].bcs[ibc, iedgelo][4] = iedgehi + 1
                    mesh3d.elem[iel + 1].bcs[ibc, 5][0] = "E"
                    mesh3d.elem[iel + 1].bcs[ibc, 5][3] = iel + 1 + 2
                    mesh3d.elem[iel + 1].bcs[ibc, 5][4] = iedgehi + 1
                    mesh3d.elem[iel + 2].bcs[ibc, iedgehi][0] = "E"
                    mesh3d.elem[iel + 2].bcs[ibc, iedgehi][3] = iel + 1 + 1
                    mesh3d.elem[iel + 2].bcs[ibc, iedgehi][4] = 6
                    mesh3d.elem[iel + 3].bcs[ibc, iedgehi][0] = "E"
                    mesh3d.elem[iel + 3].bcs[ibc, iedgehi][3] = iel + 1 + 4
                    mesh3d.elem[iel + 3].bcs[ibc, iedgehi][4] = 5
                    mesh3d.elem[iel + 4].bcs[ibc, 4][0] = "E"
                    mesh3d.elem[iel + 4].bcs[ibc, 4][3] = iel + 1 + 3
                    mesh3d.elem[iel + 4].bcs[ibc, 4][4] = iedgehi + 1
                    mesh3d.elem[iel + 4].bcs[ibc, iedgelo][0] = "E"
                    mesh3d.elem[iel + 4].bcs[ibc, iedgelo][3] = iel + 1 + 5
                    mesh3d.elem[iel + 4].bcs[ibc, iedgelo][4] = iedgehi + 1
                    mesh3d.elem[iel + 5].bcs[ibc, iedgehi][0] = "E"
                    mesh3d.elem[iel + 5].bcs[ibc, iedgehi][3] = iel + 1 + 4
                    mesh3d.elem[iel + 5].bcs[ibc, iedgehi][4] = iedgelo + 1

                # update the conditions for side faces.
                for iface in iedgelat:
                    bc = mesh.elem[i].bcs[ibc, iface][0]
                    if bc == "P" or bc == "E":
                        for iell in range(6):
                            connected_i = mesh.elem[i].bcs[ibc, iface][3] - 1
                            connected_iel = 6 * (connected_i + nel2d * (k // 4)) + iell
                            mesh3d.elem[iel + iell].bcs[ibc, iface][3] = (
                                connected_iel + 1
                            )
                            mesh3d.elem[iel + iell].bcs[ibc, iface][4] = mesh.elem[
                                i
                            ].bcs[ibc, iface][4]

    # now fix the end boundary conditions
    # face 5 is at zmin and face 6 is at zmax (with Nek indexing, corresponding to 4 and 5 in Python)
    for i in range(0, 6 * nel2d, 6):
        for ibc in range(nbc):
            i1 = i + nel3d - 6 * nel2d + 5
            mesh3d.elem[i].bcs[ibc, 4][0] = bc1[ibc]
            mesh3d.elem[i].bcs[ibc, 4][1] = i + 1
            mesh3d.elem[i].bcs[ibc, 4][2] = 5
            mesh3d.elem[i + 1].bcs[ibc, 4][0] = bc1[ibc]
            mesh3d.elem[i + 1].bcs[ibc, 4][1] = i + 1 + 1
            mesh3d.elem[i + 1].bcs[ibc, 4][2] = 5
            mesh3d.elem[i1].bcs[ibc, 5][0] = bc2[ibc]
            mesh3d.elem[i1].bcs[ibc, 5][1] = i1 + 1
            mesh3d.elem[i1].bcs[ibc, 5][2] = 6
            mesh3d.elem[i1 - 1].bcs[ibc, 5][0] = bc2[ibc]
            mesh3d.elem[i1 - 1].bcs[ibc, 5][1] = i1 - 1 + 1
            mesh3d.elem[i1 - 1].bcs[ibc, 5][2] = 6

            mesh3d.elem[i].bcs[ibc, 4][3] = 0.0
            mesh3d.elem[i].bcs[ibc, 4][4] = 0.0
            mesh3d.elem[i + 1].bcs[ibc, 4][3] = 0.0
            mesh3d.elem[i + 1].bcs[ibc, 4][4] = 0.0
            mesh3d.elem[i1].bcs[ibc, 5][3] = 0.0
            mesh3d.elem[i1].bcs[ibc, 5][4] = 0.0
            mesh3d.elem[i1 - 1].bcs[ibc, 5][3] = 0.0
            mesh3d.elem[i1 - 1].bcs[ibc, 5][4] = 0.0

            # fix the matching faces for the periodic conditions
            if bc1[ibc] == "P":
                mesh3d.elem[i].bcs[ibc, 4][3] = i1 + 1
                mesh3d.elem[i].bcs[ibc, 4][4] = 6
                mesh3d.elem[i + 1].bcs[ibc, 4][3] = i1 - 1 + 1
                mesh3d.elem[i + 1].bcs[ibc, 4][4] = 6
            if bc2[ibc] == "P":
                mesh3d.elem[i1].bcs[ibc, 5][3] = i + 1
                mesh3d.elem[i1].bcs[ibc, 5][4] = 5
                mesh3d.elem[i1 - 1].bcs[ibc, 5][3] = i + 1 + 1
                mesh3d.elem[i1 - 1].bcs[ibc, 5][4] = 5

    # FIND THE CURVED ELEMENTS
    ncurv = 0
    for el in mesh3d.elem:
        for iedge in range(12):
            if el.ccurv[iedge] != "":
                ncurv = ncurv + 1
    mesh3d.ncurv = ncurv

    # return the extruded mesh
    return mesh3d


# =================================================================================


def edge_mid(el, iedge):
    """Finds the coordinates of the midsize-node of edge iedge of element el (in other words, if the curvature were type 'm', the values of el.curv[iedge][:3]):

    Parameters
    ----------
    el : :class:`pymech.core.HexaData`
         element of mesh (usually, el=mesh.elem[i])
    iedge : int
            index of edge
    """

    # correct if ccurv=='m', otherwise, works as allocation
    midpoint = el.curv[iedge][:3]

    if el.ccurv[iedge] != "m":
        if iedge == 0:
            pos1 = el.pos[:, 0, 0, 0]
            pos2 = el.pos[:, 0, 0, 1]
        elif iedge == 1:
            pos1 = el.pos[:, 0, 0, 1]
            pos2 = el.pos[:, 0, 1, 1]
        elif iedge == 2:
            pos1 = el.pos[:, 0, 1, 1]
            pos2 = el.pos[:, 0, 1, 0]
        elif iedge == 3:
            pos1 = el.pos[:, 0, 1, 0]
            pos2 = el.pos[:, 0, 0, 0]
        elif iedge == 4:
            pos1 = el.pos[:, 1, 0, 0]
            pos2 = el.pos[:, 1, 0, 1]
        elif iedge == 5:
            pos1 = el.pos[:, 1, 0, 1]
            pos2 = el.pos[:, 1, 1, 1]
        elif iedge == 6:
            pos1 = el.pos[:, 1, 1, 1]
            pos2 = el.pos[:, 1, 1, 0]
        elif iedge == 7:
            pos1 = el.pos[:, 1, 1, 0]
            pos2 = el.pos[:, 1, 0, 0]
        elif iedge == 8:
            pos1 = el.pos[:, 0, 0, 0]
            pos2 = el.pos[:, 1, 0, 0]
        elif iedge == 9:
            pos1 = el.pos[:, 0, 0, 1]
            pos2 = el.pos[:, 1, 0, 1]
        elif iedge == 10:
            pos1 = el.pos[:, 0, 1, 1]
            pos2 = el.pos[:, 1, 1, 1]
        elif iedge == 11:
            pos1 = el.pos[:, 0, 1, 0]
            pos2 = el.pos[:, 1, 1, 0]

        if el.ccurv[iedge] == "":
            midpoint = (pos1 + pos2) / 2.0
        elif el.ccurv[iedge] == "C":
            # Curvature 'C' only needs x and y. Works for 2d and extruded meshes.
            if iedge > 7:
                # For iedge=8-11: will give a different value to what Nek considers (Nek ignores it).
                logger.warn(
                    "Calculating midpoint differently from Nek5000. Nek ignores it for edges 9-12."
                )
            radius = el.curv[iedge][0]
            dmid = (
                abs(radius)
                - (
                    radius**2
                    - (pos2[0] - pos1[0]) ** 2 / 4.0
                    - (pos2[1] - pos1[1]) ** 2 / 4.0
                )
                ** 0.5
            )
            midpoint[0] = (pos2[0] + pos1[0]) / 2.0 + dmid / (
                (pos2[0] - pos1[0]) ** 2 + (pos2[1] - pos1[1]) ** 2
            ) ** (0.5) * radius / abs(radius) * (pos2[1] - pos1[1])
            midpoint[1] = (pos2[1] + pos1[1]) / 2.0 - dmid / (
                (pos2[0] - pos1[0]) ** 2 + (pos2[1] - pos1[1]) ** 2
            ) ** (0.5) * radius / abs(radius) * (pos2[0] - pos1[0])
            midpoint[2] = (pos2[2] + pos1[2]) / 2.0
        elif el.ccurv[iedge] == "s":
            # It doesn't check if sphere is consistent with pos1 and pos2. Just assumes it is.
            radius = el.curv[iedge][4]
            center = el.curv[iedge][:3]
            dist = (pos2 + pos1) / 2.0 - center
            midpoint = (
                center
                + dist * radius / (dist[0] ** 2 + dist[1] ** 2 + dist[2] ** 2) ** 0.5
            )

    # return the coordinate of midsize node
    return midpoint


# ==============================================================================
def edge_circle(el, iedge, midpoint):
    """Finds the radius of curvature and circle center based on the midsize-node of edge iedge of element el:

    Parameters
    ----------
    el : :class:`pymech.core.HexaData`
         element of mesh (usually, el=mesh.elem[i])
    iedge : int
            index of edge
    midpoint : float
               list of coordinates of midsize-node (in other words, if the curvature were type 'm', the values of el.curv[iedge][:3])
    """

    if iedge == 0:
        pos1 = el.pos[:, 0, 0, 0]
        pos2 = el.pos[:, 0, 0, 1]
    elif iedge == 1:
        pos1 = el.pos[:, 0, 0, 1]
        pos2 = el.pos[:, 0, 1, 1]
    elif iedge == 2:
        pos1 = el.pos[:, 0, 1, 1]
        pos2 = el.pos[:, 0, 1, 0]
    elif iedge == 3:
        pos1 = el.pos[:, 0, 1, 0]
        pos2 = el.pos[:, 0, 0, 0]
    elif iedge == 4:
        pos1 = el.pos[:, 1, 0, 0]
        pos2 = el.pos[:, 1, 0, 1]
    elif iedge == 5:
        pos1 = el.pos[:, 1, 0, 1]
        pos2 = el.pos[:, 1, 1, 1]
    elif iedge == 6:
        pos1 = el.pos[:, 1, 1, 1]
        pos2 = el.pos[:, 1, 1, 0]
    elif iedge == 7:
        pos1 = el.pos[:, 1, 1, 0]
        pos2 = el.pos[:, 1, 0, 0]
    elif iedge == 8:
        pos1 = el.pos[:, 0, 0, 0]
        pos2 = el.pos[:, 1, 0, 0]
    elif iedge == 9:
        pos1 = el.pos[:, 0, 0, 1]
        pos2 = el.pos[:, 1, 0, 1]
    elif iedge == 10:
        pos1 = el.pos[:, 0, 1, 1]
        pos2 = el.pos[:, 1, 1, 1]
    elif iedge == 11:
        pos1 = el.pos[:, 0, 1, 0]
        pos2 = el.pos[:, 1, 1, 0]

    side1 = midpoint - pos1
    side2 = pos2 - midpoint
    side3 = pos1 - pos2

    d1 = (side1[0] ** 2 + side1[1] ** 2 + side1[2] ** 2) ** 0.5
    d2 = (side2[0] ** 2 + side2[1] ** 2 + side2[2] ** 2) ** 0.5
    d3 = (side3[0] ** 2 + side3[1] ** 2 + side3[2] ** 2) ** 0.5
    sper = (d1 + d2 + d3) / 2.0
    area = (sper * (sper - d1) * (sper - d2) * (sper - d3)) ** 0.5

    if area > 0.0001 * d1 * d2:
        radius = d1 * d2 * d3 / (4 * area)
        alpha1 = d2**2 * (d1**2 + d3**2 - d2**2) / 2.0
        alpha2 = d3**2 * (d2**2 + d1**2 - d3**2) / 2.0
        alpha3 = d1**2 * (d3**2 + d2**2 - d1**2) / 2.0
        center = (alpha1 * pos1 + alpha2 * midpoint + alpha3 * pos2) / (8.0 * area**2)
        if (
            (side1[0] - side3[0]) * (side2[1] - side1[1])
            - (side1[0] - side2[0]) * (side3[1] - side1[1])
        ) < 0.0:
            # if curvature == 'C', the radius is negative for clockwise triangles
            # works only for 2d/extruded mesh - do not know how to interpret it in 3d (should work for edges edges 0-7 of extruded meshes, unknown behaviour for edges 8-11)
            radius = -radius
    else:
        # radius is too big compared to edge. For pratical purposes, no curvature: radius=0.0
        radius = 0.0
        center = [0.0, 0.0, 0.0]

    curv = el.curv[iedge][:4]
    curv[0] = radius
    curv[1:4] = center

    # return radius and center of curvature (it is not ready to be used: For 'C', define curv[1:4]=0.0; For 's', define el.curv[4]=abs(curv[0]) and el.curv[:3]=curv[1:4])
    return curv


# =================================================================================


def delete_internal_bcs(mesh):
    """Deletes the internal boundary conditions 'E' in a Nek5000 mesh.
    Those are present in .rea files but do not need to be valid, and are completely absent from .re2 files.
    Returns the number of deleted conditions.

    Parameters
    ----------
    mesh : :class:`pymech.core.HexaData`
           The mesh to modify in-place.
    """

    ndelete = 0
    for ibc in range(mesh.nbc):
        for el in mesh.elem:
            for iface in range(6):
                bc = el.bcs[ibc, iface]
                if bc[0] == "E":
                    ndelete = ndelete + 1
                    bc[0] = ""
                    for i in range(1, 8):
                        bc[i] = 0
    return ndelete


# =================================================================================


def generate_internal_bcs(mesh, tol=1e-3):
    """Generates internal boundary conditions 'E' in a Nek5000 mesh based on geometric proximity of faces.
    This will loop over all pairs of faces and connect them if their centres are closer than some specified relative tolerance.
    This function returns the number of internal connections found.

    Parameters
    ----------
    mesh : :class:`pymech.core.HexaData`
           The mesh to modify in-place

    tol: float
         the tolerance within which the centres of two faces are considered the same,
         relative to the smallest edge of the elements
    """

    # First generate a length scale for each element, equal to the smallest edge of that element.
    scales = np.zeros((mesh.nel,))
    for iel, el in enumerate(mesh.elem):
        scales[iel] = el.smallest_edge()

        # check if there is a zero length edge; in this case the mesh is invalid and there is no point continuing.
        if scales[iel] <= 0.0:
            raise ValueError(f"Detected an edge with zero length in element {iel}!")

    # generate lookup tables for face centers and the faces that are already connected
    nface = 2 * mesh.ndim
    fc = np.empty((mesh.nel, nface, 3))
    lc = np.full((mesh.nel, nface), False, dtype=bool)
    for iel in range(mesh.nel):
        el = mesh.elem[iel]
        for iface in range(nface):
            fc[iel, iface, :] = el.face_center(iface)
            if el.bcs[0, iface][0] != "":
                lc[iel, iface] = True  # mark connected faces as True

    # Now that we have the scales, we can compare the location of the faces for each pair of elements and connect them if they are close
    nconnect = 0  # number of connections made
    for iel1 in range(mesh.nel):
        el1 = mesh.elem[iel1]
        lmin1 = scales[iel1]
        for iface1 in range(nface):
            if not lc[iel1, iface1]:  # skip if face already connected
                xf0, yf0, zf0 = fc[iel1, iface1, :]
                find_iel1_iface = False  # for break criterion
                for iel2 in range(iel1 + 1, mesh.nel):
                    el2 = mesh.elem[iel2]
                    lmin2 = scales[iel2]
                    max_d = tol * min(lmin1, lmin2)
                    for iface2 in range(nface):
                        if not lc[iel2, iface2]:  # skip if face already connected
                            xf1, yf1, zf1 = fc[iel2, iface2, :]
                            dist = np.sqrt(
                                (xf1 - xf0) ** 2 + (yf1 - yf0) ** 2 + (zf1 - zf0) ** 2
                            )
                            if dist <= max_d:
                                for ibc in range(mesh.nbc):
                                    # increment counter for diagnostics
                                    nconnect = nconnect + 1
                                    # write the connectivity information in both directions
                                    el1.bcs[ibc, iface1][0] = "E"
                                    el1.bcs[ibc, iface1][1] = iel1 + 1
                                    el1.bcs[ibc, iface1][2] = iface1 + 1
                                    el1.bcs[ibc, iface1][3] = iel2 + 1
                                    el1.bcs[ibc, iface1][4] = iface2 + 1
                                    el2.bcs[ibc, iface2][0] = "E"
                                    el2.bcs[ibc, iface2][1] = iel2 + 1
                                    el2.bcs[ibc, iface2][2] = iface2 + 1
                                    el2.bcs[ibc, iface2][3] = iel1 + 1
                                    el2.bcs[ibc, iface2][4] = iface1 + 1
                                # update logical arrays
                                lc[iel1, iface1] = True
                                lc[iel2, iface2] = True
                                # each face of iel1 can only be connected to a single
                                # face of another element, break loops once connected
                                find_iel1_iface = True
                                break
                    if find_iel1_iface:
                        break  # this is the long loop iel2 over all elements!
    return nconnect


# =================================================================================


def keep_elements(mesh: HexaData, elems, external_bc=""):
    """
    Reduce the mesh to a subset of its elements

    Parameters
    ----------
    mesh : :class:`pymech.core.HexaData`
            mesh to modify
    elems: int array
            list of element numbers (zero indexed) to keep
    external_bc: string
            if an element was connected to another which is removed from the mesh, replace that condition with external_bc
    """

    # sort element numbers
    elems = np.array(elems)
    elems.sort()

    # list of booleans representing whether each element is kept
    kept = [False for i in range(mesh.nel)]

    # list of offsets to apply to each element index
    offsets = np.zeros((mesh.nel,))
    last_iel = -1
    current_offset = 0
    for iel in elems:
        if iel >= mesh.nel:
            raise ValueError(f"invalid element number {iel} for nel = {mesh.nel}")
        current_offset += iel - last_iel - 1
        offsets[iel] = current_offset
        last_iel = iel
        kept[iel] = True

    # apply negative offsets to internal & periodic boundary conditions
    for iel, ibc, iface in product(elems, range(mesh.nbc), range(2 * mesh.ndim)):
        offset = offsets[iel]
        bc = mesh.elem[iel].bcs[ibc, iface][0]
        mesh.elem[iel].bcs[ibc, iface][1] = iel - offset + 1
        if bc == "E" or bc == "P":
            connected_iel = int(mesh.elem[iel].bcs[ibc, iface][3]) - 1
            if kept[connected_iel]:
                # update the index of the connected element
                mesh.elem[iel].bcs[ibc, iface][3] -= offsets[connected_iel]
            else:
                # delete the boundary condition or replace it with an external one
                mesh.elem[iel].bcs[ibc, iface][0] = external_bc
                mesh.elem[iel].bcs[ibc, iface][3] = 0
                mesh.elem[iel].bcs[ibc, iface][4] = 0

    # delete unused elements from element list
    new_elems = []
    for iel in elems:
        new_elems.append(mesh.elem[iel])
    mesh.elem = new_elems
    mesh.nel = len(mesh.elem)
    mesh.update_ncurv()


# =================================================================================


def exponential_refinement_parameter(l0, ltot, n, tol=1e-14):
    """
    In a 1D exponential mesh spacing where n elements of lengths l0, alpha*l0, ..., alpha^(n-1)*l0
    add up to a total length ltot, return the parameter alpha given l0, ltot, and n.

    Parameters
    ----------
    l0 : float
        the initial mesh spacing (size of the first element)
    ltot : float
        total length of the mesh
    n : integer
        number of elements
    tol : float
        error tolerance on alpha

    Returns
    -------
    alpha : float
        the mesh refinement parameter
    """

    # interval around 1 where to approximate the functions (see below)
    eps_approx_f = 2.0e-6
    eps_approx_fp = 2.0e-6
    # maximum number of iterations; really this should converge in less than 100.
    maxiter = 10000

    # the total length is ltot = (1 - alpha^n)/(1 - alpha) * l0,
    # so we want to find alpha such that (1 - alpha^n)/(1 - alpha) = ltot/l0.
    # It's better posed numerically if we solve instead for log((1 - alpha^n)/(1 - alpha)) = log(ltot/l0),
    def err(x):
        # We can further approximate the function around x=1 to get rid of the 0/0 appearing in the expression.
        # In practice it seems to be worth using the linearised form for |alpha-1| < 2e-6 or so for n=1000 and l/l0=10000
        if abs(x - 1) < eps_approx_f:
            # this would be first order accurate
            # return log(n * l0 / ltot) + 0.5 * (n - 1) * (x - 1)
            # this is third order
            return log(n * l0 / ltot) + log(
                1
                + (n - 1) / 2 * (x - 1)
                + (n - 1) * (n - 2) / 6 * (x - 1) ** 2
                + (n - 1) * (n - 2) * (n - 3) / 24 * (x - 1) ** 3
            )
        else:
            return log((x**n - 1) / (x - 1) * l0 / ltot)

    def err_prime(x):
        # We can further approximate the derivative of the error function around x=1
        # to get rid of the 0^2/0^2 appearing in the expression.
        # In practice it seems to be worth using the linearised form for |alpha-1| < 2e-6 or so for n=1000 and l/l0=10000
        if abs(x - 1) < eps_approx_fp:
            return 0.5 * (n - 1) + (n - 1) * (n - 5) / 12 * (x - 1)
        else:
            # this would be first order accurate
            # return n * x**(n - 1) / (x**n - 1) - 1 / (x - 1)
            # this is second order
            return (
                (n - 1) / 2
                + (n - 1) * (n - 2) / 3 * (x - 1)
                + (n - 1) * (n - 2) * (n - 3) / 8 * (x - 1) ** 2
            ) / (1 + (n - 1) / 2 * (x - 1) + (n - 1) * (n - 2) / 6 * (x - 1) ** 2)

    # Solve this using the Newton method.
    # initial guess: we can choose alpha = 1 now that the function is well behaved around that point.
    alpha = 1
    residual = 1
    iter = 0
    while residual > tol and iter < maxiter:
        fx = err(alpha)
        fpx = err_prime(alpha)
        alpha1 = alpha - fx / fpx
        residual = abs(alpha1 - alpha)
        alpha = alpha1
        logger.debug(f"{iter} {alpha}, res: {residual}")
        iter += 1
    return alpha


# =================================================================================


def rotate_2d(mesh, x0, y0, theta):
    """
    Rotate a mesh around an axis aligned with z passing through (x0, y0, 0) by an angle theta.

    Parameters
    ----------
    mesh : :class:`pymech.core.HexaData`
        the mesh to modify in place
    x0 : float
        x-coordinate of the center of rotation
    y0 : float
        y-coordinate of the center of rotation
    theta : rotation angle
    """

    cost = cos(theta)
    sint = sin(theta)
    for el in mesh.elem:
        x = el.pos[0, ...].copy()
        y = el.pos[1, ...].copy()
        el.pos[0, ...] = cost * x - sint * y
        el.pos[1, ...] = sint * x + cost * y
        # rotate 'm' curvature
        for edge in range(12):
            if el.ccurv[edge] == "m":
                x = el.curv[edge, 0]
                y = el.curv[edge, 1]
                el.curv[edge, 0] = cost * x - sint * y
                el.curv[edge, 1] = sint * x + cost * y


# =================================================================================


def gen_box(
    nx: int,
    ny: int,
    xmin: float,
    xmax: float,
    ymin: float,
    ymax: float,
    var=[2, 2, 1, 0, 0],
    bcs_xmin=None,
    bcs_xmax=None,
    bcs_ymin=None,
    bcs_ymax=None,
    internal_bcs=True,
):
    """
    generates a rectangular box of nx×ny elements in [xmin, xmax]×[ymin, ymax].
    Boundary conditions can optionally be scpecified for each side.

    Parameters
    ----------
    nx : float
          number of spectral elements in the x direction
    ny : float
          number of spectral elements in the y direction
    xmin, xmax, ymin, ymax : float
          bounds pof the box in the x and y directions
    var : int list
          optional, the list of fields to include in the mesh (velocity, pressure, temperature, passive scalars).
          Defaults to only velocity and pressure.
    bcs_xmin, bcs_xmax, bcs_ymin, bcs_ymax : str lists
          optional, lists of boundary conditions for each of the sides of the box.
          Each argument must be either `None` or a list of boundary conditions, one for each field present in the mesh except pressure. If `None`, it defaults to periodic boundary conditions for all fields.
    internal_bcs : bool
          Optional, specifies whether to build the connectivity boundary conditions between internal elements. Default is `True`.
    """

    lr1 = [2, 2, 1]  # the mesh is 2D so lz1 = 1
    ndim = 2
    nel = nx * ny
    # count the number of boundary conditions required
    nbc = 0
    if var[1] > 0:
        nbc += 1
    if var[3] > 0:
        nbc += 1
    nbc += var[4]

    # set default periodic boundary conditions
    if bcs_xmin is None:
        bcs_xmin = ["P"] * nbc
    if bcs_xmax is None:
        bcs_xmax = ["P"] * nbc
    if bcs_ymin is None:
        bcs_ymin = ["P"] * nbc
    if bcs_ymax is None:
        bcs_ymax = ["P"] * nbc

    box = HexaData(ndim, nel, lr1, var, nbc=nbc)
    # box dimensions
    lx = xmax - xmin
    if lx <= 0:
        raise ValueError(f"xmax must be greater than xmin, but xmax - xmin = {lx:9e}")
    ly = ymax - ymin
    if ly <= 0:
        raise ValueError(f"ymax must be greater than ymin, but ymax - ymin = {ly:9e}")

    # indexing in each block, 0-indexed
    def elnum(i, j, ni, nj):
        return i + ni * j

    for i, j in product(range(nx), range(ny)):
        # coordinates of corners
        x0 = xmin + lx * i / nx
        x1 = xmin + lx * (i + 1) / nx
        y0 = ymin + ly * j / ny
        y1 = ymin + ly * (j + 1) / ny
        # assign the coordinates
        el = box.elem[elnum(i, j, nx, ny)]
        el.pos[0, 0, :, 0] = x0
        el.pos[0, 0, :, 1] = x1
        el.pos[1, 0, 0, :] = y0
        el.pos[1, 0, 1, :] = y1
    # connectivity
    if internal_bcs:
        for ibc, i, j in product(range(nbc), range(nx), range(ny)):
            el = box.elem[elnum(i, j, nx, ny)]
            # bottom face
            if j != 0:
                el.bcs[ibc, 0][0] = "E"
                el.bcs[ibc, 0][1] = elnum(i, j, nx, ny) + 1
                el.bcs[ibc, 0][2] = 1
                el.bcs[ibc, 0][3] = elnum(i, j - 1, nx, ny) + 1
                el.bcs[ibc, 0][4] = 3
            # right face
            if i != nx - 1:
                el.bcs[ibc, 1][0] = "E"
                el.bcs[ibc, 1][1] = elnum(i, j, nx, ny) + 1
                el.bcs[ibc, 1][2] = 2
                el.bcs[ibc, 1][3] = elnum(i + 1, j, nx, ny) + 1
                el.bcs[ibc, 1][4] = 4
            # top face
            if j != ny - 1:
                el.bcs[ibc, 2][0] = "E"
                el.bcs[ibc, 2][1] = elnum(i, j, nx, ny) + 1
                el.bcs[ibc, 2][2] = 3
                el.bcs[ibc, 2][3] = elnum(i, j + 1, nx, ny) + 1
                el.bcs[ibc, 2][4] = 1
            # left face
            if i != 0:
                el.bcs[ibc, 3][0] = "E"
                el.bcs[ibc, 3][1] = elnum(i, j, nx, ny) + 1
                el.bcs[ibc, 3][2] = 4
                el.bcs[ibc, 3][3] = elnum(i - 1, j, nx, ny) + 1
                el.bcs[ibc, 3][4] = 2

    # Apply boundary conditions
    for ibc in range(nbc):
        for i in range(nx):
            # bottom face
            j = 0
            el = box.elem[elnum(i, j, nx, ny)]
            bc = bcs_ymin[ibc]
            el.bcs[ibc, 0][0] = bc
            el.bcs[ibc, 0][1] = elnum(i, j, nx, ny) + 1
            el.bcs[ibc, 0][2] = 1
            if bc == "P":
                el.bcs[ibc, 0][3] = elnum(i, ny - 1, nx, ny) + 1
                el.bcs[ibc, 0][4] = 3
            # top face
            j = ny - 1
            el = box.elem[elnum(i, j, nx, ny)]
            bc = bcs_ymax[ibc]
            el.bcs[ibc, 2][0] = bc
            el.bcs[ibc, 2][1] = elnum(i, j, nx, ny) + 1
            el.bcs[ibc, 2][2] = 3
            if bc == "P":
                el.bcs[ibc, 2][3] = elnum(i, 0, nx, ny) + 1
                el.bcs[ibc, 2][4] = 1
        for j in range(ny):
            # left face
            i = 0
            el = box.elem[elnum(i, j, nx, ny)]
            bc = bcs_xmin[ibc]
            el.bcs[ibc, 3][0] = bc
            el.bcs[ibc, 3][1] = elnum(i, j, nx, ny) + 1
            el.bcs[ibc, 3][2] = 4
            if bc == "P":
                el.bcs[ibc, 3][3] = elnum(nx - 1, j, nx, ny) + 1
                el.bcs[ibc, 3][4] = 2
            # right face
            i = nx - 1
            el = box.elem[elnum(i, j, nx, ny)]
            bc = bcs_xmax[ibc]
            el.bcs[ibc, 1][0] = bc
            el.bcs[ibc, 1][1] = elnum(i, j, nx, ny) + 1
            el.bcs[ibc, 1][2] = 2
            if bc == "P":
                el.bcs[ibc, 1][3] = elnum(0, j, nx, ny) + 1
                el.bcs[ibc, 1][4] = 4

    return box


# =================================================================================


def gen_circle(
    r: float,
    s: float,
    ns: int,
    no: int,
    curvature_fun=None,
    bl_fun=None,
    var=[2, 2, 1, 0, 0],
    bc=["W"],
    internal_bcs=True,
):
    """
    Generates a 2D circular mesh with a square at the center surrounded by an O mesh.

    Parameters
    ----------
    r : float
        radius of the mesh
    s : float
        relative length of the diagonal of the centre square mesh to the circle diameter
    ns : int
        number of elements in the side of the square and in a quarter of the circumference of the circle
    no : int
        number of elements in the O mesh part in the radial direction
    internal_bcs : bool
        if True, builds the internal connectivity information of the elements
    curvature_fun : [0, 1] -> [0, 1] function or None
        Function defining the evolution of the relative curvature of the edges between the
        straight side of the square and the circular edge.
        Defaults to x -> sin(pi/2 x) if None. A constant 1 means concentric circles,
        a constant zero means only straight edges inside the domain, but ideally you want
        a continuous function with f(0) = 0 and f(1) = 1.
    bl_fun : [0, 1] -> [0, 1] function or None
        Function defining the evolution of the grid location between the edge of
        the square at 0 and the edge of the circle at 1.
        Defaults to an exponential grid spacing with a uniform spacing around the square if None.
    var : integer list
        Number of geometry, velocity, pressure, temperature and scalar variables to define in the mesh
    bc : str list
        boundary conditions to use for each field
    internal_bcs : bool
        whether to build internal connectivity information
    """

    # Generate five sub-meshes: one square and four quarter O meshes.
    # Because of the symmetry, generate only one quarter and rotate it.

    # dimension constants for mesh / elements generation
    lr1 = [2, 2, 1]  # the mesh is 2D so lz1 = 1
    ndim = 2
    nbc = len(bc)

    # default curvature function
    if curvature_fun is None:

        def curvature_fun(x):
            return sin(0.5 * pi * x) ** 2
            # return sqrt(1 - (1 - x)**2)  # this one also works

    # default curvature and BL functions
    if bl_fun is None:
        # element size in the square
        square_spacing = s * r * sqrt(2) / ns
        # total length between the side of the square to the edge of the largest square fitting in the circle
        width = 0.5 * r * (1 - s) * sqrt(2)
        # we want the size of the first element to be `square_spacing`, but it's going to be stretched by the curvature function,
        # so we need to find the length before stretching that gives the right length after stretching

        def err_stretching(len_new, len_ref):
            eta = curvature_fun(len_new / width)
            total_width = r * (1 - 0.5 * sqrt(2) * s)
            l1 = len_new / width * ((1 - eta) * width + eta * total_width)
            return l1 - len_ref

        # iterate using the secant method
        l0 = 1
        l1 = 0
        f0 = err_stretching(l0, square_spacing)
        eps = 1e-12  # this value shouldn't matter very much
        max_iters = 10000
        iters = 0
        while abs(f0) > eps and iters < max_iters:
            f1 = err_stretching(l1, square_spacing)
            # Update l0 and l1
            l0, l1 = l1, l1 - (l1 - l0) / (f1 - f0) * f1
            f0 = f1
            iters += 1

        # geometric spacing progression factor using that first element size
        alpha = exponential_refinement_parameter(l1, width, no)

        def bl_fun(x):
            eps = 1e-12
            # (l0 + alpha*l0 + ... + alpha^j*l0)/ltot = (l0 * (alpha^j - 1)/(alpha - 1))/ltot
            # with j = x*no
            if abs(alpha - 1) < eps:
                return x
            else:
                return l1 / width * (alpha ** (x * no) - 1) / (alpha - 1)

    # indexing in each block, 0-indexed
    def elnum(i, j, ni, nj):
        return i + ni * j

    # external boundary conditions
    def apply_bcs(mesh, ni, nj, bc0, bc1, bc2, bc3):
        for ibc in range(nbc):
            for i in range(ni):
                # bottom face
                mesh.elem[elnum(i, 0, ni, nj)].bcs[ibc, 0][0] = bc0[ibc]
                # top face
                mesh.elem[elnum(i, nj - 1, ni, nj)].bcs[ibc, 2][0] = bc2[ibc]
            for j in range(nj):
                # right face
                mesh.elem[elnum(ni - 1, j, ni, nj)].bcs[ibc, 1][0] = bc1[ibc]
                # left face
                mesh.elem[elnum(0, j, ni, nj)].bcs[ibc, 3][0] = bc3[ibc]

    # mesh connectivity for a structured grid, ignoring boundaries
    def build_connectivity(mesh, ni, nj):
        for i, j in product(range(ni), range(nj)):
            el = mesh.elem[elnum(i, j, ni, nj)]
            for ibc in range(nbc):
                # bottom face
                if j != 0:
                    el.bcs[ibc, 0][0] = "E"
                    el.bcs[ibc, 0][1] = elnum(i, j, ni, nj) + 1
                    el.bcs[ibc, 0][2] = 1
                    el.bcs[ibc, 0][3] = elnum(i, j - 1, ni, nj) + 1
                    el.bcs[ibc, 0][4] = 3
                # right face
                if i != ni - 1:
                    el.bcs[ibc, 1][0] = "E"
                    el.bcs[ibc, 1][1] = elnum(i, j, ni, nj) + 1
                    el.bcs[ibc, 1][2] = 2
                    el.bcs[ibc, 1][3] = elnum(i + 1, j, ni, nj) + 1
                    el.bcs[ibc, 1][4] = 4
                # top face
                if j != nj - 1:
                    el.bcs[ibc, 2][0] = "E"
                    el.bcs[ibc, 2][1] = elnum(i, j, ni, nj) + 1
                    el.bcs[ibc, 2][2] = 3
                    el.bcs[ibc, 2][3] = elnum(i, j + 1, ni, nj) + 1
                    el.bcs[ibc, 2][4] = 1
                # left face
                if i != 0:
                    el.bcs[ibc, 3][0] = "E"
                    el.bcs[ibc, 3][1] = elnum(i, j, ni, nj) + 1
                    el.bcs[ibc, 3][2] = 4
                    el.bcs[ibc, 3][3] = elnum(i - 1, j, ni, nj) + 1
                    el.bcs[ibc, 3][4] = 2

    # Box 1: square
    nel_square = ns * ns
    box_square = HexaData(ndim, nel_square, lr1, var, nbc=nbc)
    # the characteristic distance is the diagonal
    c = s * r * sqrt(2)
    for i, j in product(range(ns), range(ns)):
        # coordinates of corners
        x0 = c * (i / ns - 0.5)
        x1 = c * ((i + 1) / ns - 0.5)
        y0 = c * (j / ns - 0.5)
        y1 = c * ((j + 1) / ns - 0.5)
        # make a new element
        el = box_square.elem[elnum(i, j, ns, ns)]
        el.pos[0, 0, :, 0] = x0
        el.pos[0, 0, :, 1] = x1
        el.pos[1, 0, 0, :] = y0
        el.pos[1, 0, 1, :] = y1
    # connectivity
    if internal_bcs:
        build_connectivity(box_square, ns, ns)
    # boundary conditions: dummy BCs to signal that the faces should be glued
    connectivity_bc = ["con"] * nbc
    apply_bcs(
        box_square,
        ns,
        ns,
        connectivity_bc,
        connectivity_bc,
        connectivity_bc,
        connectivity_bc,
    )

    # Box 2: quarter-O
    nel_o = no * ns
    box_o = HexaData(ndim, nel_o, lr1, var, nbc=nbc)
    for i, j in product(range(no), range(ns)):
        # angular positions, between -pi/4 and pi/4
        alpha0 = 0.5 * pi * (j / ns - 0.5)
        alpha1 = 0.5 * pi * ((j + 1) / ns - 0.5)
        # angular positions along the square
        beta0 = atan(2 * j / ns - 1)
        beta1 = atan(2 * (j + 1) / ns - 1)
        # distance from centre of circle to square edge
        l_sa0 = 0.5 * c / cos(alpha0)
        l_sa1 = 0.5 * c / cos(alpha1)
        l_sb0 = 0.5 * c / cos(beta0)
        l_sb1 = 0.5 * c / cos(beta1)
        # distance from square edge to circle edge
        l_c0 = r - l_sa0
        l_c1 = r - l_sa1
        # distance from square edge to edge of largest square fitting in the circle
        l_h0 = 0.5 * sqrt(2) * r / cos(beta0) - l_sb0
        l_h1 = 0.5 * sqrt(2) * r / cos(beta1) - l_sb1
        # progression along the radius
        eta0 = bl_fun(i / no)
        eta1 = bl_fun((i + 1) / no)
        # interpolation between l_h and l_c using curvature function
        gamma0 = curvature_fun(eta0)
        gamma1 = curvature_fun(eta1)
        d_c00 = l_sa0 + eta0 * l_c0
        d_c01 = l_sa0 + eta1 * l_c0
        d_c10 = l_sa1 + eta0 * l_c1
        d_c11 = l_sa1 + eta1 * l_c1
        d_h00 = l_sb0 + eta0 * l_h0
        d_h01 = l_sb0 + eta1 * l_h0
        d_h10 = l_sb1 + eta0 * l_h1
        d_h11 = l_sb1 + eta1 * l_h1
        # angles from center
        theta00 = gamma0 * alpha0 + (1 - gamma0) * beta0
        theta01 = gamma1 * alpha0 + (1 - gamma1) * beta0
        theta10 = gamma0 * alpha1 + (1 - gamma0) * beta1
        theta11 = gamma1 * alpha1 + (1 - gamma1) * beta1
        # distances of corners to center
        d00 = gamma0 * d_c00 + (1 - gamma0) * d_h00
        d01 = gamma1 * d_c01 + (1 - gamma1) * d_h01
        d10 = gamma0 * d_c10 + (1 - gamma0) * d_h10
        d11 = gamma1 * d_c11 + (1 - gamma1) * d_h11
        # coordinates of corners
        x00 = d00 * cos(theta00)
        x01 = d01 * cos(theta01)
        x10 = d10 * cos(theta10)
        x11 = d11 * cos(theta11)
        y00 = d00 * sin(theta00)
        y01 = d01 * sin(theta01)
        y10 = d10 * sin(theta10)
        y11 = d11 * sin(theta11)
        # make a new element
        el = box_o.elem[elnum(i, j, no, ns)]
        el.pos[0, 0, 0, 0] = x00
        el.pos[0, 0, 0, 1] = x01
        el.pos[0, 0, 1, 0] = x10
        el.pos[0, 0, 1, 1] = x11
        el.pos[1, 0, 0, 0] = y00
        el.pos[1, 0, 0, 1] = y01
        el.pos[1, 0, 1, 0] = y10
        el.pos[1, 0, 1, 1] = y11
    # connectivity
    if internal_bcs:
        build_connectivity(box_o, no, ns)
    # boundary conditions: dummy BCs on the faces to be connected, external BC on the right face
    apply_bcs(box_o, no, ns, connectivity_bc, bc, connectivity_bc, connectivity_bc)
    # add circular curvature for external faces
    for j in range(ns):
        el = box_o.elem[elnum(no - 1, j, no, ns)]
        edge = 1  # right edge
        el.ccurv[edge] = "C"
        el.curv[edge, 0] = r

    # copy the O box, rotate it to the other sides and merge it into the mesh
    # right box
    box_square.merge(box_o, ignore_all_bcs=not internal_bcs)
    # top box
    rotate_2d(box_o, 0, 0, 0.5 * pi)
    box_square.merge(box_o, ignore_all_bcs=not internal_bcs)
    # left box
    rotate_2d(box_o, 0, 0, 0.5 * pi)
    box_square.merge(box_o, ignore_all_bcs=not internal_bcs)
    # bottom box
    rotate_2d(box_o, 0, 0, 0.5 * pi)
    box_square.merge(box_o, ignore_all_bcs=not internal_bcs)

    return box_square


# =================================================================================


def map2D(
    mesh: HexaData,
    transformation,
    curvature=True,
    boundary_curvature=True,
):
    """
    Applies a coordinate transformation to a 2D mesh, returning the transformed mesh.
    The faces are curved to second-order accuracy by setting the midpoints according to the transformation.


    Parameters
    ----------
    mesh : :class:`pymech.core.HexaData`
        the mesh to transform, will not be modified.
    transformation : (float, float) -> (float float) function
        the coordinate transformation to apply. It must be a valid right-handed transformation (the determinant of its Jacobian must be positive on the mesh domain).
    curvature : bool
        specifies whether to apply curvature to all faces. The faces that are already curved in the original mesh will always be curved. Curvature can still be applied on the boundaries when this is `False` if `boundary_curvature` is set to `True`. True by default.
    boundary_curvature: bool
        specifies whether to apply curvature to boundary faces by setting midpoints. `True` by default.

    Returns
    -------
    mapped_mesh : :class:`pymech.core.HexaData`
        the transformed mesh
    """

    mapped_mesh = copy.deepcopy(mesh)
    for el, mapped_el in zip(mesh.elem, mapped_mesh.elem):
        # map the vertices of the elements. This is the easy part.
        # we might want to use vectorisation here but the input function might not support it
        for ix, iy in product(range(2), range(2)):
            x = el.pos[0, 0, ix, iy]
            y = el.pos[1, 0, ix, iy]
            mapped_x, mapped_y = transformation(x, y)
            mapped_el.pos[0, 0, ix, iy] = mapped_x
            mapped_el.pos[1, 0, ix, iy] = mapped_y

        # Now, fix the curvature.
        # Several scenarios are possible:
        # - the edge is curved with a midpoint. In this case, we update the midpoint.
        # - the edge is curved with a circle arc ("C"). In this case, we transform this into a midpoint curvature and update it, because there is no guarantee that the image of a circle by the transformation is a circle.
        # - the edge is a boundary edge and boundary_curvature=True or any edge and curvature=True. then we generate a new midpoint at the image of the centre of the original edge.
        # otherwise, leave the edge without curvature.
        for iedge in range(4):
            if (
                curvature
                or el.ccurv[iedge] != ""
                or (boundary_curvature and el.bcs[0][iedge][0] not in ["", "E"])
            ):
                mapped_el.ccurv[iedge] = "m"
                xm, ym, _ = edge_mid(el, iedge)
                mapped_xm, mapped_ym = transformation(xm, ym)
                mapped_el.curv[iedge, 0:2] = mapped_xm, mapped_ym
    return mapped_mesh
