import numpy as np
import pymech.exadata as exdat
import copy
from pymech.log import logger
from itertools import product


# ==============================================================================
def extrude(mesh: exdat, z, bc1='P', bc2='P', internal_bcs=True):
    """Extrudes a 2D mesh into a 3D one

    Parameters
    ----------
    mesh : exadata
           2D mesh structure to extrude
    zmin : float
           min value of the z coordinate to extrude to
    zmax : float
           max value of the z coordinate to extrude to
    z :  float 1d array
           z coordinates at which to extrude the mesh
    bc1, bc2 : str
           the boundary conditions to use at each end
    internal_bcs : bool
           if True, build mesh connectivity using internal 'E' boundary conditions
           (note that those are not used by Nek5000 and will not be written to binary .re2 files).
    """

    if mesh.ndim != 2:
        logger.critical('The mesh to extrude must be 2D')
        return -1
    if mesh.lr1 != [2, 2, 1]:
        logger.critical('Only mesh structures can be extruded (lr1 = [2, 2, 1])')
        return -2
    if mesh.var[0] < 2:
        logger.critical('The mesh to extrude must contain (x, y) geometry')
        return -3
    if (bc1 == 'P' and bc2 != 'P') or (bc1 != 'P' and bc2 == 'P'):
        logger.critical('Inconsistent boundary conditions: one end is \'P\' but the other isn\'t')
        return -4

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
        mesh3d.elem = mesh3d.elem + copy.deepcopy(mesh.elem)

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
                mesh3d.elem[iel].curv[icurv + 4] = mesh3d.elem[iel].curv[icurv]    # curvature params
                if curv_type == 'm':
                    # in this case the midpoint is given. (x, y) is correct but z should be set to the proper value.
                    mesh3d.elem[iel].curv[icurv][2] = z1
                    mesh3d.elem[iel].curv[icurv + 4][2] = z2

    # fix the internal boundary conditions (even though it's probably useless)
    # the end boundary conditions will be overwritten later with the proper ones
    if internal_bcs:
        for (iel, el) in enumerate(mesh3d.elem):
            for ibc in range(nbc):
                el.bcs[ibc, 4][0] = 'E'
                el.bcs[ibc, 4][1] = iel + 1
                el.bcs[ibc, 4][2] = 5
                el.bcs[ibc, 4][3] = iel - nel2d + 1
                el.bcs[ibc, 4][4] = 6
                el.bcs[ibc, 5][0] = 'E'
                el.bcs[ibc, 5][1] = iel + 1
                el.bcs[ibc, 5][2] = 6
                el.bcs[ibc, 5][3] = iel + nel2d + 1
                el.bcs[ibc, 5][4] = 5
                # update the conditions for side faces
                for iface in range(4):
                    el.bcs[ibc, iface][1] = iel + 1
                    if el.bcs[ibc, iface][0] == 'E':
                        # el.bcs[ibc, 0][1] ought to contain iel+1 once the mesh is valid
                        # but for now it should be off by a factor of nel2d because it is a copy of an element in the first slice
                        offset = iel - el.bcs[ibc, iface][1] + 1
                        el.bcs[ibc, iface][3] = el.bcs[ibc, iface][3] + offset

    # now fix the end boundary conditions
    # face 5 is at zmin and face 6 is at zmax (with Nek indexing, corresponding to 4 and 5 in Python)
    for i in range(nel2d):
        for ibc in range(nbc):
            i1 = i + (nz - 1) * nel2d  # index of the face on the zmax side
            mesh3d.elem[i].bcs[ibc, 4][0] = bc1
            mesh3d.elem[i].bcs[ibc, 4][1] = i + 1
            mesh3d.elem[i].bcs[ibc, 4][2] = 5
            mesh3d.elem[i1].bcs[ibc, 5][0] = bc2
            mesh3d.elem[i1].bcs[ibc, 5][1] = i1 + 1
            mesh3d.elem[i1].bcs[ibc, 5][2] = 6
            # fix the matching faces for the periodic conditions
            if bc1 == 'P':
                mesh3d.elem[i].bcs[ibc, 4][3] = i1 + 1
                mesh3d.elem[i].bcs[ibc, 4][4] = 6
            if bc2 == 'P':
                mesh3d.elem[i1].bcs[ibc, 5][3] = i + 1
                mesh3d.elem[i1].bcs[ibc, 5][4] = 5

    # return the extruded mesh
    return mesh3d


# ==============================================================================
def extrude_refine(mesh2D, z, bc1='P', bc2='P', fun=None, funpar=None, imesh_high=0):
    """Extrudes a 2D mesh into a 3D one, following the pattern:
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
    mesh2D : exadata
           2D mesh structure to extrude
    z : float array
        list of z values of the  most refined zones of the extruded mesh
    bc : str
         the boundary condition to use at the ends
    fun: function
         list of functions that define the splitting lines for different discretization meshes (default: empty, in which case the simple extrusion function `extrude` is called instead)
    funpar: list
          list of parameters for functions that define the splitting lines for different discretization meshes (default: None, equivalent to an array of zeroes)
    imesh_high : int
                 index of fun that defines the mesh with higher discretization. Example: 0, is the most internal mesh; 1 is the second most internal mesh, etc (default: the most internal mesh, imesh_high=0)
    """

    # Consistency checks: Initial grid
    if mesh2D.ndim != 2:
        logger.critical('The mesh to extrude must be 2D')
        return -1
    if mesh2D.lr1 != [2, 2, 1]:
        logger.critical('Only mesh structures can be extruded (lr1 = [2, 2, 1])')
        return -2
    if mesh2D.var[0] < 2:
        logger.critical('The mesh to extrude must contain (x, y) geometry')
        return -3
    # Consistency checks: Periodic boundary condition
    if (bc1 == 'P' and bc2 != 'P') or (bc1 != 'P' and bc2 == 'P'):
        logger.critical('Inconsistent boundary conditions: one end is \'P\' but the other isn\'t')
        return -4

    # Consistency checks: Functions that define the splitting lines
    nsplit = len(fun)
    if funpar is not None and len(funpar) != nsplit:
        logger.critical(f'The length of funpar ({len(funpar)}) must match the length of par ({nsplit})!')
        return -5

    # number of elements in the z direction
    nz = len(z) - 1

    # Consistency checks: if nz is divided by 4 (or 8, or 16, etc)
    if (nz % 2**abs(imesh_high + 1) != 0) or (nz % 2**abs(nsplit - imesh_high + 1) != 0):
        logger.critical(f'Inconsistent elements to extrude: the number of elements ({nz}) must be a multiple of {max([2**abs(imesh_high + 1), 2**abs(nsplit - imesh_high + 1)])}')
        return -10

    # If fun is not defined, there is no splitting to be done. Call simple extrusion and end routine
    if fun is None:
        logger.info('Splitting function not defined. Calling simple extrusion routine.')
        mesh3D = extrude(mesh2D, z, bc1, bc2)
        return mesh3D

    mesh2D_ext = copy.deepcopy(mesh2D)
    # get rid of internal boundary conditions because they will be hard to get correct and they are not used by Nek5000 anyway
    delete_internal_bcs(mesh2D_ext)

    meshes2D = []  # list of 2D meshes
    meshes3D = []  # list of 3D meshes

    # Splitting 2D meshes
    # sorts element of the 2D mesh into layers: most refined, first split, second most refined,
    # second split, etc, until the least refined stored into mesh2D_ext.
    for k in range(nsplit):
        meshes2D.append(copy.deepcopy(mesh2D_ext))
        meshes2D.append(copy.deepcopy(mesh2D_ext))

        iel_int = 0
        iel_mid = 0
        iel_ext = 0

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
            if max(rvec) <= 0.0:  # the element belongs to the internal (more refined) mesh
                meshes2D[2 * k].elem[iel_int] = copy.deepcopy(mesh2D_ext.elem[iel])
                iel_int += 1
            elif min(rvec) > 0.0:  # the element belongs to the external (less refined) mesh
                mesh2D_ext.elem[iel_ext] = copy.deepcopy(mesh2D_ext.elem[iel])
                iel_ext += 1
            else:  # the element belongs to the intermediate mesh and will be split
                meshes2D[2 * k + 1].elem[iel_mid] = copy.deepcopy(mesh2D_ext.elem[iel])
                iel_mid += 1

        meshes2D[2 * k].nel = iel_int
        meshes2D[2 * k + 1].nel = iel_mid
        mesh2D_ext.nel = iel_ext

        logger.debug(f'Mesh2D {2 * k} elements: {meshes2D[2 * k].nel}')
        logger.debug(f'Mesh2D {2 * k + 1} elements: {meshes2D[2 * k + 1].nel}')

        meshes2D[2 * k].elem = meshes2D[2 * k].elem[:iel_int]
        meshes2D[2 * k + 1].elem = meshes2D[2 * k + 1].elem[:iel_mid]
        mesh2D_ext.elem = mesh2D_ext.elem[:iel_ext]

    # End of splitting, remaining is the last mesh: Mesh_ext
    logger.debug(f'Mesh2Dext elements: {mesh2D_ext.nel}')

    # Update curvature metadata for all sub-meshes
    mesh2D_ext.update_ncurv()
    for mesh_part in meshes2D:
        mesh_part.update_ncurv()

    # Extruding meshes
    logger.info('Extruding meshes')
    for k in range(nsplit):
        # divide number of elements by 2**(coarsening level)
        n_local = nz // 2**abs(k - imesh_high)
        # select z coordinates for coarsened elements
        z_local = z[::int(2**abs(k - imesh_high))]

        if k < imesh_high:
            def fun_local(xpos, ypos, rlim):
                return -fun[k](xpos, ypos, rlim)
            n_mid = 2 * n_local
            z_mid = z[::int(2**abs(k - imesh_high + 1))]
        else:
            fun_local = fun[k]
            n_mid = n_local
            z_mid = z_local

        if n_mid % 4 != 0:
            logger.critical(f'Inconsistent elements to extrude: n ({n_mid}) is not a multiple of 4.')
            return -11

        meshes3D.append(extrude(meshes2D[2 * k], z_local, bc1=bc1, bc2=bc2, internal_bcs=False))
        meshes3D.append(extrude_mid(meshes2D[2 * k + 1], z_mid, bc1, bc2, fun_local, funpar[k]))

        logger.debug(f'Mesh3D {2 * k} elements: {meshes3D[2 * k].nel}')
        logger.debug(f'Mesh3D {2 * k + 1} elements: {meshes3D[2 * k + 1].nel}')

    n_local = nz // 2**abs(nsplit - imesh_high)
    z_local = z[::int(2**abs(nsplit - imesh_high))]
    mesh3D_ext = extrude(mesh2D_ext, z_local, bc1=bc1, bc2=bc2, internal_bcs=False)

    logger.debug(f'Mesh3Dext elements: {mesh3D_ext.nel}')

    # Merging meshes
    logger.info('Merging meshes')
    mesh3D = mesh3D_ext
    for mesh_part in meshes3D:
        mesh3D.merge(mesh_part, ignore_all_bcs=True)
    logger.info(f'Merging done. Total elements: {mesh3D.nel}')

    # update curve metadata
    mesh3D.update_ncurv()

    # return the extruded mesh
    return mesh3D


# =================================================================================

def extrude_mid(mesh, z, bc1, bc2, fun, funpar=0.0):
    """Extrudes the mid elments of the 2D mesh into a 3D one. Following the pattern:
     _____ _____ _____ _____
    |1   /|\   4|    /|\    |
    |__ / | \ __|__ / | \ __| (fun (with parameter funpar) should change change sign in the mid element)
    |0 |2 |3 | 5|  |  |  |  | (half of the mid elements are also divided in 2 in (x, y)-plane)
    |__|__|__|__|__|__|__|__| (numbers in the figure indicate the indices (iel + 0; iel + 1; etc))

    Parameters
    ----------
    mesh : exadata
           2D mesh structure to extrude
    z : float
        list of z values of the nodes of the elements of the extruded mesh in the high discretization region (len(z)-1 must be divide by 4)
    bc : str
         the boundary condition to use at the ends
    fun : function
          function that define the splitting lines for different discretization meshes
    funpar : not defined, depends on the function
             parameter for functions that define the splitting lines for different discretization meshes (default: zero, can be used for when funpar is not needed inside fun)

    Suggestion: see function extrude_split to understand how to call extrude_mid
    """

    # Consistency checks: Initial grid
    if mesh.ndim != 2:
        logger.critical('The mesh to extrude must be 2D')
        return -1
    if mesh.lr1 != [2, 2, 1]:
        logger.critical('Only mesh structures can be extruded (lr1 = [2, 2, 1])')
        return -2
    if mesh.var[0] < 2:
        logger.critical('The mesh to extrude must contain (x, y) geometry')
        return -3
    if (bc1 == 'P' and bc2 != 'P') or (bc1 != 'P' and bc2 == 'P'):
        logger.critical('Inconsistent boundary conditions: one end is \'P\' but the other isn\'t')
        return -4

    nz = len(z) - 1
    z1 = np.zeros((nz, 1))
    z2 = np.zeros((nz, 1))
    z1 = z[0:nz]
    z2 = z[1:nz + 1]

    if nz % 4 != 0:
        logger.critical('Inconsistent elements to extrude: nz must be divided by 4')
        return -5

    # copy the structure and make it 3D
    mesh3d = copy.deepcopy(mesh)
    mesh3d.lr1 = [2, 2, 2]
    mesh3d.var = [3, 0, 0, 0, 0]  # remove anything that isn't geometry for now
    nel2d = mesh.nel
    nel3d = nel2d * 6 * (nz // 4)  # every mid-extrusion of a 2d-element creates 6 3d-elements, while the usual extrusion creates 4 in the high discretized grid and 2 in the low discretized grid
    nbc = mesh.nbc
    mesh3d.nel = nel3d
    mesh3d.ndim = 3
    # The curved sides will also be extruded, one on each side of each element along nz
    mesh3d.ncurv = 2 * nz * mesh.ncurv

    # add extra copies of all elements
    for k in range(6 * (nz // 4) - 1):
        mesh3d.elem = mesh3d.elem + copy.deepcopy(mesh.elem)

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
            index_lo = np.zeros((2, 2), dtype=int)  # index of low points
            index_hi = np.zeros((2, 2), dtype=int)  # index of high points
            iindex_lo = 0
            iindex_hi = 0
            iedgelat = np.zeros((2), dtype=int)
            iedgeconlat = np.zeros((2), dtype=int)
            for ii in range(2):
                for jj in range(2):
                    xvec[jj, ii] = mesh.elem[i].pos[0, 0, jj, ii]
                    yvec[jj, ii] = mesh.elem[i].pos[1, 0, jj, ii]
                    rvec[jj, ii] = fun(xvec[jj, ii], yvec[jj, ii], funpar)
                    if rvec[jj, ii] <= 0.0:
                        if (iindex_lo > 1):
                            logger.critical('Mid element not consistent. Criteria must divide elements with 2 points on each side.')
                            return -11
                        index_lo[iindex_lo, :] = [jj, ii]
                        iindex_lo += 1
                    else:
                        if (iindex_hi > 1):
                            logger.critical('Mid element not consistent. Criteria must divide elements with 2 points on each side.')
                            return -11
                        index_hi[iindex_hi, :] = [jj, ii]
                        iindex_hi += 1
            if (iindex_lo != 2) or (iindex_hi != 2):
                logger.critical('Mid element not consistent. Criteria must divide elements with 2 points on each side.')
                return -11

            # find the indices of edges, for curvature and boundary condition
            # iedgehi is the index of the edge of element iel + 0 that is the intersection between iel + 0 and iel + 1 (high edge). iedgelo is the index of the edge of element iel + 1 that is the intersection (low edge).
            # iedgelat are the indices of the lateral (splitted) edges. iedgeconlat are the indices of the edges (edge in z-direction) of elements iel + 2 and iel + 3 that connect to the respective lateral edges.
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
            poslo = copy.deepcopy(mesh.elem[i].pos[:, :1, index_lo[:, 0], index_lo[:, 1]])
            poshi = copy.deepcopy(mesh.elem[i].pos[:, :1, index_hi[:, 0], index_hi[:, 1]])
            # mid position is influenced by curvature
            posmid = 0.5 * (mesh.elem[i].pos[:, :1, index_lo[:, 0], index_lo[:, 1]] + mesh.elem[i].pos[:, :1, index_hi[:, 0], index_hi[:, 1]])
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
                if mesh.elem[i].ccurv[icurv] == 's':
                    logger.warn(f'Curvature s on element {i + 1}. Not consistent with extrusion, changing to C')
                    mesh.elem[i].ccurv[icurv] == 'C'
                    mesh.elem[i].ccurv[icurv][0] = mesh.elem[i].ccurv[icurv][4]
                    mesh.elem[i].ccurv[icurv][1:4] = 0.0
                elif (mesh.elem[i].ccurv[icurv] != '') and (mesh.elem[i].ccurv[icurv] != 'm') and (mesh.elem[i].ccurv[icurv] != 'C'):
                    logger.warn(f'Curvature unknown on element {i + 1}')

            # extend curvature and correct it if necessary
            # calculate coordinates of midsize-node for every curvature type, except both empty (even if both are 'C'). Then find radius, if applicable. 'm' takes precendence over 'C'.
            # if mesh.elem[i].ccurv[iedgehi] != mesh.elem[i].ccurv[iedgelo]:
            if mesh.elem[i].ccurv[iedgehi] != '' or mesh.elem[i].ccurv[iedgelo] != '':
                midpointhi = edge_mid(mesh.elem[i], iedgehi)
                midpointlo = edge_mid(mesh.elem[i], iedgelo)
                midpointmid = 0.5 * (posmid[:, 0, 0] + posmid[:, 0, 1]) + 0.5 * (midpointhi - 0.5 * (poshi[:, 0, 0] + poshi[:, 0, 1]) + midpointlo - 0.5 * (poslo[:, 0, 0] + poslo[:, 0, 1]))
                if mesh.elem[i].ccurv[iedgehi] == 'm' or mesh.elem[i].ccurv[iedgelo] == 'm':
                    mesh3d.elem[iel].ccurv[iedgehi] = 'm'
                    mesh3d.elem[iel + 1].ccurv[iedgelo] = 'm'
                    mesh3d.elem[iel].curv[iedgehi][:3] = midpointmid
                    mesh3d.elem[iel + 1].curv[iedgelo][:3] = midpointmid
                elif mesh.elem[i].ccurv[iedgehi] == 'C' or mesh.elem[i].ccurv[iedgelo] == 'C':
                    midpointmid[2] = z1[k]
                    curvmid = edge_circle(mesh3d.elem[iel], iedgehi, midpointmid)
                    if curvmid[0] == 0.0:
                        mesh3d.elem[iel].ccurv[iedgehi] = ''
                        mesh3d.elem[iel + 1].ccurv[iedgelo] = ''
                        mesh3d.elem[iel].curv[iedgehi][:4] = 0.0
                        mesh3d.elem[iel + 1].curv[iedgelo][:4] = 0.0
                    else:
                        curvmid[1:4] = 0.0
                        mesh3d.elem[iel].ccurv[iedgehi] = 'C'
                        mesh3d.elem[iel + 1].ccurv[iedgelo] = 'C'
                        mesh3d.elem[iel].curv[iedgehi][:4] = curvmid
                        mesh3d.elem[iel + 1].curv[iedgelo][:4] = -curvmid
                else:
                    logger.warn(f'Splitted element curvature unknown on element {i + 1} of mid mesh. Removing curvature.')
                    # For cases not implemented, remove curvature
                    mesh3d.elem[iel].ccurv[iedgehi] = ''
                    mesh3d.elem[iel + 1].ccurv[iedgelo] = ''
                    mesh3d.elem[iel].curv[iedgehi] = 0.0 * (mesh.elem[i].curv[iedgehi])
                    mesh3d.elem[iel + 1].curv[iedgelo] = copy.deepcopy(mesh3d.elem[iel].curv[iedgehi])

            for ilat in range(2):
                # Fixing curvature of edges divided in half. For curv_type == 'C', it is not a true extrusion - 'diagonal edges' not consistent with 'C'.
                if mesh.elem[i].ccurv[iedgelat[ilat]] == 'm':
                    # coordinates of midsize-node is approximated. Considering an edge aligned with the x-axis, the position would be (ym_new = 3/4*ym_old = 1/2*ym_old(mean y-position) + 1/4*ym_old(distance)) at xm_new=(x2-x1)/4. This comes from parabolic curve (however, it is not exactly midsize)
                    dmid = ((poshi[0, 0, ilat] - poslo[0, 0, ilat]) * (poslo[1, 0, ilat] - posmid[1, 0, ilat]) - (poslo[0, 0, ilat] - posmid[0, 0, ilat]) * (poshi[1, 0, ilat] - poslo[1, 0, ilat])) / ((poshi[0, 0, ilat] - poslo[0, 0, ilat])**2 + (poshi[1, 0, ilat] - poslo[1, 0, ilat])**2)
                    mesh3d.elem[iel].curv[iedgelat[ilat]][0] = 0.5 * (poslo[0, 0, ilat] + posmid[0, 0, ilat]) + dmid / 4.0 * (poshi[1, 0, ilat] - poslo[1, 0, ilat])
                    mesh3d.elem[iel + 1].curv[iedgelat[ilat]][0] = 0.5 * (posmid[0, 0, ilat] + poshi[0, 0, ilat]) + dmid / 4.0 * (poshi[1, 0, ilat] - poslo[1, 0, ilat])
                    mesh3d.elem[iel].curv[iedgelat[ilat]][1] = 0.5 * (poslo[1, 0, ilat] + posmid[1, 0, ilat]) - dmid / 4.0 * (poshi[0, 0, ilat] - poslo[0, 0, ilat])
                    mesh3d.elem[iel + 1].curv[iedgelat[ilat]][1] = 0.5 * (posmid[1, 0, ilat] + poshi[1, 0, ilat]) - dmid / 4.0 * (poshi[0, 0, ilat] - poslo[0, 0, ilat])
                elif mesh.elem[i].ccurv[iedgelat[ilat]] == 'C':
                    # if the lateral edge has curvature 'C', the diagonal edges connected to it (iedgeconlat of elements iel + 2 and iel + 3) would have curvature 'C', which are inconsistent with edges 8-12 inside Nek5000 (because these are supposed to be edges in z-direction). Changing to curvature 'm' for elements iel + 1 (and iel + 2, iel + 3, iel + 2, iel + 4)
                    midpointlathi = edge_mid(mesh3d.elem[iel + 1], iedgelat[ilat])
                    mesh3d.elem[iel + 1].curv[iedgelat[ilat]][:3] = midpointlathi
                    mesh3d.elem[iel + 1].ccurv[iedgelat[ilat]] = 'm'

            for icurv in range(4):
                # a 2D element has 4 edges that can be curved, numbered 0-3;
                # the extruded 3D element can have four more (on the other side), numbered 4-7
                for iell in range(6):
                    mesh3d.elem[iel + iell].ccurv[icurv + 4] = mesh3d.elem[iel + iell].ccurv[icurv]
                    mesh3d.elem[iel + iell].curv[icurv + 4] = mesh3d.elem[iel + iell].curv[icurv]    # curvature params

            mesh3d.elem[iel + 2].curv[0:4] = copy.deepcopy(mesh3d.elem[iel].curv[0:4])
            mesh3d.elem[iel + 3].curv[4:8] = copy.deepcopy(mesh3d.elem[iel].curv[0:4])
            mesh3d.elem[iel + 5].curv = copy.deepcopy(mesh3d.elem[iel].curv)
            mesh3d.elem[iel + 4].curv = copy.deepcopy(mesh3d.elem[iel + 1].curv)
            mesh3d.elem[iel + 2].ccurv[0:4] = copy.deepcopy(mesh3d.elem[iel].ccurv[0:4])
            mesh3d.elem[iel + 3].ccurv[4:8] = copy.deepcopy(mesh3d.elem[iel].ccurv[0:4])
            mesh3d.elem[iel + 5].ccurv = copy.deepcopy(mesh3d.elem[iel].ccurv)
            mesh3d.elem[iel + 4].ccurv = copy.deepcopy(mesh3d.elem[iel + 1].ccurv)

            for icurv in range(4):
                # z should be set to the proper value.
                curv_type = mesh3d.elem[iel].ccurv[icurv]
                if curv_type == 'm':
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
                if curv_type == 'm':
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
                        mesh3d.elem[iel + 1].curv[icurv + 4][izcurv] = 0.5 * (z1[k + 1] + z2[k + 1])
                        mesh3d.elem[iel + 4].curv[icurv][izcurv] = 0.5 * (z2[k + 1] + z2[k + 2])

            # Fixing the curvature of 3d-edges in z-direction that connects to lateral edges in trapezoidal elements (all other edges in z-direction - indices 8 to 11 - do not have curvature)
            for ilat in range(2):
                mesh3d.elem[iel + 2].curv[iedgeconlat[ilat]] = copy.deepcopy(mesh3d.elem[iel + 1].curv[iedgelat[ilat] + 4])
                mesh3d.elem[iel + 2].ccurv[iedgeconlat[ilat]] = copy.deepcopy(mesh3d.elem[iel + 1].ccurv[iedgelat[ilat] + 4])
                mesh3d.elem[iel + 3].curv[iedgeconlat[ilat]] = copy.deepcopy(mesh3d.elem[iel + 4].curv[iedgelat[ilat]])
                mesh3d.elem[iel + 3].ccurv[iedgeconlat[ilat]] = copy.deepcopy(mesh3d.elem[iel + 4].ccurv[iedgelat[ilat]])

    # fix the internal boundary conditions (even though it's probably useless)
    # the end boundary conditions will be overwritten later with the proper ones
            for ibc in range(nbc):
                # set the conditions in faces normal to z
                for iell in range(6):
                    mesh3d.elem[iel + iell].bcs[ibc, 4][0] = 'E'
                    mesh3d.elem[iel + iell].bcs[ibc, 4][1] = iel + iell + 1
                    mesh3d.elem[iel + iell].bcs[ibc, 4][2] = 5
                    mesh3d.elem[iel + iell].bcs[ibc, 4][4] = 6
                    mesh3d.elem[iel + iell].bcs[ibc, 5][0] = 'E'
                    mesh3d.elem[iel + iell].bcs[ibc, 5][1] = iel + iell + 1
                    mesh3d.elem[iel + iell].bcs[ibc, 5][2] = 6
                    mesh3d.elem[iel + iell].bcs[ibc, 5][4] = 5

                mesh3d.elem[iel].bcs[ibc, 4][3] = iel + 1 - nel2d
                mesh3d.elem[iel].bcs[ibc, 5][3] = iel + 1 + 2
                mesh3d.elem[iel + 1].bcs[ibc, 4][3] = iel + 1 - nel2d
                mesh3d.elem[iel + 1].bcs[ibc, 5][3] = iel + 1 + 2
                mesh3d.elem[iel + 2].bcs[ibc, 4][3] = iel + 1
                mesh3d.elem[iel + 2].bcs[ibc, 5][3] = iel + 1 + 3
                mesh3d.elem[iel + 3].bcs[ibc, 4][3] = iel + 1 + 2
                mesh3d.elem[iel + 3].bcs[ibc, 5][3] = iel + 1 + 5
                mesh3d.elem[iel + 4].bcs[ibc, 4][3] = iel + 1 + 3
                mesh3d.elem[iel + 4].bcs[ibc, 5][3] = iel + 1 + nel2d
                mesh3d.elem[iel + 5].bcs[ibc, 4][3] = iel + 1 + 3
                mesh3d.elem[iel + 5].bcs[ibc, 5][3] = iel + 1 + nel2d
                # update the conditions for side faces. (FIXME : Not corrected for merging - need to know numbering of other simply-extruded meshes (it is not really necessary, internal bc are not used))
                for iface in range(4):
                    for iell in range(6):
                        mesh3d.elem[iel + iell].bcs[ibc, iface][1] = iel + iell + 1
                        mesh3d.elem[iel + iell].bcs[ibc, iface][2] = iface + 1
                    if mesh3d.elem[iel].bcs[ibc, iface][0] == 'E':
                        # el.bcs[ibc, 0][1] ought to contain iel + 1 once the mesh is valid
                        # but for now it should be off by a factor of nel2d because it is a copy of an element in the first slice
                        ielneigh = 6 * (mesh3d.elem[iel].bcs[ibc, iface][3] - 1 + nel2d * (k // 4))
                        for iell in range(6):
                            mesh3d.elem[iel + iell].bcs[ibc, iface][3] = ielneigh + 1 + iell

                # Correct internal bc for mid faces of elements.
                mesh3d.elem[iel].bcs[ibc, iedgehi][0] = 'E'
                mesh3d.elem[iel].bcs[ibc, iedgehi][3] = iel + 1 + 1
                mesh3d.elem[iel].bcs[ibc, iedgehi][4] = iedgelo + 1
                mesh3d.elem[iel + 1].bcs[ibc, iedgelo][0] = 'E'
                mesh3d.elem[iel + 1].bcs[ibc, iedgelo][3] = iel + 1
                mesh3d.elem[iel + 1].bcs[ibc, iedgelo][4] = iedgehi + 1
                mesh3d.elem[iel + 1].bcs[ibc, 5][0] = 'E'
                mesh3d.elem[iel + 1].bcs[ibc, 5][3] = iel + 1 + 2
                mesh3d.elem[iel + 1].bcs[ibc, 5][4] = iedgehi + 1
                mesh3d.elem[iel + 2].bcs[ibc, iedgehi][0] = 'E'
                mesh3d.elem[iel + 2].bcs[ibc, iedgehi][3] = iel + 1 + 1
                mesh3d.elem[iel + 2].bcs[ibc, iedgehi][4] = 6
                mesh3d.elem[iel + 3].bcs[ibc, iedgehi][0] = 'E'
                mesh3d.elem[iel + 3].bcs[ibc, iedgehi][3] = iel + 1 + 4
                mesh3d.elem[iel + 3].bcs[ibc, iedgehi][4] = 5
                mesh3d.elem[iel + 4].bcs[ibc, 4][0] = 'E'
                mesh3d.elem[iel + 4].bcs[ibc, 4][3] = iel + 1 + 3
                mesh3d.elem[iel + 4].bcs[ibc, 4][4] = iedgehi + 1
                mesh3d.elem[iel + 4].bcs[ibc, iedgelo][0] = 'E'
                mesh3d.elem[iel + 4].bcs[ibc, iedgelo][3] = iel + 1 + 5
                mesh3d.elem[iel + 4].bcs[ibc, iedgelo][4] = iedgehi + 1
                mesh3d.elem[iel + 5].bcs[ibc, iedgehi][0] = 'E'
                mesh3d.elem[iel + 5].bcs[ibc, iedgehi][3] = iel + 1 + 4
                mesh3d.elem[iel + 5].bcs[ibc, iedgehi][4] = iedgelo + 1

    # now fix the end boundary conditions
    # face 5 is at zmin and face 6 is at zmax (with Nek indexing, corresponding to 4 and 5 in Python)
    for i in range(0, 6 * nel2d, 6):
        for ibc in range(nbc):
            i1 = i + nel3d - 6 * nel2d + 5
            mesh3d.elem[i].bcs[ibc, 4][0] = bc1
            mesh3d.elem[i].bcs[ibc, 4][1] = i + 1
            mesh3d.elem[i].bcs[ibc, 4][2] = 5
            mesh3d.elem[i + 1].bcs[ibc, 4][0] = bc1
            mesh3d.elem[i + 1].bcs[ibc, 4][1] = i + 1 + 1
            mesh3d.elem[i + 1].bcs[ibc, 4][2] = 5
            mesh3d.elem[i1].bcs[ibc, 5][0] = bc2
            mesh3d.elem[i1].bcs[ibc, 5][1] = i1 + 1
            mesh3d.elem[i1].bcs[ibc, 5][2] = 6
            mesh3d.elem[i1 - 1].bcs[ibc, 5][0] = bc2
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
            if bc1 == 'P':
                mesh3d.elem[i].bcs[ibc, 4][3] = i1 + 1
                mesh3d.elem[i].bcs[ibc, 4][4] = 6
                mesh3d.elem[i + 1].bcs[ibc, 4][3] = i1 - 1 + 1
                mesh3d.elem[i + 1].bcs[ibc, 4][4] = 6
            if bc2 == 'P':
                mesh3d.elem[i1].bcs[ibc, 5][3] = i + 1
                mesh3d.elem[i1].bcs[ibc, 5][4] = 5
                mesh3d.elem[i1 - 1].bcs[ibc, 5][3] = i + 1 + 1
                mesh3d.elem[i1 - 1].bcs[ibc, 5][4] = 5

    # Removing internal boundary conditions. 'E' boundary conditions are wrong where meshes merge, but should be right internally. (When possible: FIXME indices and delete these lines. However, it is not really necessary, internal bc are not used)
    for (iel, el) in enumerate(mesh3d.elem):
        for ibc in range(mesh3d.nbc):
            for iface in range(6):
                el.bcs[ibc, iface][1] = iel + 1
                if el.bcs[ibc, iface][0] == 'E':
                    el.bcs[ibc, iface][0] = ''
                    el.bcs[ibc, iface][1] = iel + 1
                    el.bcs[ibc, iface][2] = iface + 1
                    el.bcs[ibc, iface][3] = 0.0
                    el.bcs[ibc, iface][4] = 0.0

    # FIND THE CURVED ELEMENTS
    ncurv = 0
    for el in mesh3d.elem:
        for iedge in range(12):
            if el.ccurv[iedge] != '':
                ncurv = ncurv + 1
    mesh3d.ncurv = ncurv

    # return the extruded mesh
    return mesh3d


# =================================================================================

def edge_mid(el, iedge):
    """Finds the coordinates of the midsize-node of edge iedge of element el (in other words, if the curvature were type 'm', the values of el.curv[iedge][:3]):

    Parameters
    ----------
    el : exadata
         element of mesh (usually, el=mesh.elem[i])
    iedge : int
            index of edge
    """

    # correct if ccurv=='m', otherwise, works as allocation
    midpoint = copy.deepcopy(el.curv[iedge][:3])

    if el.ccurv[iedge] != 'm':
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

        if el.ccurv[iedge] == '':
            midpoint = (pos1 + pos2) / 2.0
        elif el.ccurv[iedge] == 'C':
            # Curvature 'C' only needs x and y. Works for 2d and extruded meshes.
            if iedge > 7:
                # For iedge=8-11: will give a different value to what Nek considers (Nek ignores it).
                logger.warn('Calculating midpoint differently from Nek5000. Nek ignores it for edges 9-12.')
            radius = el.curv[iedge][0]
            dmid = abs(radius) - (radius**2 - (pos2[0] - pos1[0])**2 / 4.0 - (pos2[1] - pos1[1])**2 / 4.0)**0.5
            midpoint[0] = (pos2[0] + pos1[0]) / 2.0 + dmid / ((pos2[0] - pos1[0])**2 + (pos2[1] - pos1[1])**2)**(0.5) * radius / abs(radius) * (pos2[1] - pos1[1])
            midpoint[1] = (pos2[1] + pos1[1]) / 2.0 - dmid / ((pos2[0] - pos1[0])**2 + (pos2[1] - pos1[1])**2)**(0.5) * radius / abs(radius) * (pos2[0] - pos1[0])
            midpoint[2] = (pos2[2] + pos1[2]) / 2.0
        elif el.ccurv[iedge] == 's':
            # It doesn't check if sphere is consistent with pos1 and pos2. Just assumes it is.
            radius = el.curv[iedge][4]
            center = el.curv[iedge][:3]
            dist = (pos2 + pos1) / 2.0 - center
            midpoint = center + dist * radius / (dist[0]**2 + dist[1]**2 + dist[2]**2)**0.5

    # return the coordinate of midsize node
    return midpoint


# ==============================================================================
def edge_circle(el, iedge, midpoint):
    """Finds the radius of curvature and circle center based on the midsize-node of edge iedge of element el:

    Parameters
    ----------
    el : exadata
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

    d1 = (side1[0]**2 + side1[1]**2 + side1[2]**2)**0.5
    d2 = (side2[0]**2 + side2[1]**2 + side2[2]**2)**0.5
    d3 = (side3[0]**2 + side3[1]**2 + side3[2]**2)**0.5
    sper = (d1 + d2 + d3) / 2.0
    area = (sper * (sper - d1) * (sper - d2) * (sper - d3))**0.5

    if area > 0.0001 * d1 * d2:
        radius = d1 * d2 * d3 / (4 * area)
        alpha1 = d2**2 * (d1**2 + d3**2 - d2**2) / 2.0
        alpha2 = d3**2 * (d2**2 + d1**2 - d3**2) / 2.0
        alpha3 = d1**2 * (d3**2 + d2**2 - d1**2) / 2.0
        center = (alpha1 * pos1 + alpha2 * midpoint + alpha3 * pos2) / (8.0 * area**2)
        if ((side1[0] - side3[0]) * (side2[1] - side1[1]) - (side1[0] - side2[0]) * (side3[1] - side1[1])) < 0.0:
            # if curvature == 'C', the radius is negative for clockwise triangles
            # works only for 2d/extruded mesh - do not know how to interpret it in 3d (should work for edges edges 0-7 of extruded meshes, unknown behaviour for edges 8-11)
            radius = -radius
    else:
        # radius is too big compared to edge. For pratical purposes, no curvature: radius=0.0
        radius = 0.0
        center = [0.0, 0.0, 0.0]

    curv = copy.deepcopy(el.curv[iedge][:4])
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
    mesh : exadata
           The mesh to modify in-place.
    """

    ndelete = 0
    for ibc in range(mesh.nbc):
        for el in mesh.elem:
            for iface in range(6):
                bc = el.bcs[ibc, iface]
                if bc[0] == 'E':
                    ndelete = ndelete + 1
                    bc[0] = ''
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
    mesh : exadata
           The mesh to modify in-place

    tol: float
         the tolerance within which the centres of two faces are considered the same,
         relative to the smallest edge of the elements
    """

    # First generate a length scale (squared) for each element, equal to the square of the smallest edge of that element.
    scales = np.zeros((mesh.nel,))
    for (iel, (el, l2)) in enumerate(zip(mesh.elem, scales)):
        # get coordinates of points
        x1, y1, z1 = el.pos[:, 0, 0, 0]
        x2, y2, z2 = el.pos[:, 0, 0, 1]
        x3, y3, z3 = el.pos[:, 0, 1, 0]
        x4, y4, z4 = el.pos[:, 0, 1, 1]
        # compute squares of edges lengths
        edges_l2 = np.zeros((12,))
        edges_l2[0] = (x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2
        edges_l2[1] = (x3 - x2)**2 + (y3 - y2)**2 + (z3 - z2)**2
        edges_l2[2] = (x4 - x3)**2 + (y4 - y3)**2 + (z4 - z3)**2
        edges_l2[3] = (x1 - x4)**2 + (y1 - y4)**2 + (z1 - z4)**2
        if mesh.ndim > 2:
            # in 3D, do the same for the upper face, and also the side edges
            x5, y5, z5 = el.pos[:, 1, 0, 0]
            x6, y6, z6 = el.pos[:, 1, 0, 1]
            x7, y7, z7 = el.pos[:, 1, 1, 0]
            x8, y8, z8 = el.pos[:, 1, 1, 1]
            edges_l2[4] = (x6 - x5)**2 + (y6 - y5)**2 + (z6 - z5)**2
            edges_l2[5] = (x7 - x6)**2 + (y7 - y6)**2 + (z7 - z6)**2
            edges_l2[6] = (x8 - x7)**2 + (y8 - y7)**2 + (z8 - z7)**2
            edges_l2[7] = (x5 - x8)**2 + (y5 - y8)**2 + (z5 - z8)**2
            edges_l2[8] = (x5 - x1)**2 + (y5 - y1)**2 + (z5 - z1)**2
            edges_l2[9] = (x6 - x2)**2 + (y6 - y2)**2 + (z6 - z2)**2
            edges_l2[10] = (x7 - x3)**2 + (y7 - y3)**2 + (z7 - z3)**2
            edges_l2[11] = (x8 - x4)**2 + (y8 - y4)**2 + (z8 - z4)**2
            l2 = edges_l2.min()
        else:
            l2 = edges_l2[:4].min()

        # check if there is a zero length edge; in this case the mesh is invalid and there is no point continuing.
        if l2 <= 0.0:
            logger.critical(f'Detected an edge with zero length in element {iel}!')
            return -1

    # Now that we have the scales, we can compare the location of the faces for each pair of elements and connect them if they are close
    nconnect = 0  # number of connections made
    for iel in range(mesh.nel):
        el = mesh.elem[iel]
        l2 = scales[iel]
        for other_iel in range(iel + 1, mesh.nel):
            other_el = mesh.elem[other_iel]
            other_l2 = scales[other_iel]
            max_d2 = tol**2 * min(l2, other_l2)
            for iface in range(2 * mesh.ndim):
                xf, yf, zf = el.face_center(iface)
                for other_iface in range(2 * mesh.ndim):
                    other_xf, other_yf, other_zf = other_el.face_center(other_iface)
                    dist2 = (other_xf - xf)**2 + (other_yf - yf)**2 + (other_zf - zf)**2
                    if dist2 <= max_d2:
                        for ibc in range(mesh.nbc):
                            # increment counter for diagnostics
                            nconnect = nconnect + 1
                            # write the connectivity information in both directions
                            el.bcs[ibc, iface][0] = 'E'
                            el.bcs[ibc, iface][1] = iel + 1
                            el.bcs[ibc, iface][2] = iface + 1
                            el.bcs[ibc, iface][3] = other_iel + 1
                            el.bcs[ibc, iface][4] = other_iface + 1
                            other_el.bcs[ibc, other_iface][0] = 'E'
                            other_el.bcs[ibc, other_iface][1] = other_iel + 1
                            other_el.bcs[ibc, other_iface][2] = other_iface + 1
                            other_el.bcs[ibc, other_iface][3] = iel + 1
                            other_el.bcs[ibc, other_iface][4] = iface + 1

    return nconnect


# =================================================================================

def offset_connectivity(mesh: exdat, offset: int, iel_min=0):
    """
    Adds a value to the index of the elements connected via internal or periodic
    boundary conditions to elements of the mesh. This is used to keep the connectivity
    valid when deleting or inserting elements in the mesh.

    Parameters
    ----------
    mesh    : exadata
           The mesh to modify
    offset  : int
           The value by which to offset the indices
    iel_min : int
           The first element (in zero-based indexing) to offset
    """

    for el, ibc, iface in product(mesh.elem, range(mesh.nbc), range(2 * mesh.dim)):
        bc = el.bcs[ibc, iface][0]
        if bc == 'E' or bc == 'P':
            if int(el.bcs[ibc, iface][3]) > iel_min:  # the connected element number is 1-indexed
                el.bcs[ibc, iface][3] += offset
