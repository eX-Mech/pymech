import numpy as np
import pymech.exadata as exdat
import copy
from pymech.log import logger


# ==============================================================================
def extrude(mesh: exdat, zmin, zmax, nz, bc1='P', bc2='P'):
    """Extrudes a 2D mesh into a 3D one

    Parameters
    ----------
    mesh : exadata
           2D mesh structure to extrude
    zmin : float
                 min value of the z coordinate to extrude to
    zmax : float
                 max value of the z coordinate to extrude to
    nz : int
         number of elements to create in the z direction
    bc : str
         the boundary condition to use at the ends
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
    mesh3d.lr1 = [2, 2, 2]
    mesh3d.var = [3, 0, 0, 0, 0]  # remove anything that isn't geometry for now
    nel2d = mesh.nel
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
            z1 = zmin + k / nz * (zmax - zmin)
            z2 = zmin + (k + 1) / nz * (zmax - zmin)
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
