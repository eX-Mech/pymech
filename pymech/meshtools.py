import numpy as np
import pymech.exadata as exdat
import copy

#==============================================================================
def extrude(mesh, zmin, zmax, nz, bc1='P', bc2='P'):
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
    nel3d = mesh.nel*nz
    nbc = mesh.nbc
    mesh3d.nel = nel3d
    mesh3d.ndim = 3

    # add extra copies of all elements
    for k in range(nz-1):
        mesh3d.elem = mesh3d.elem + copy.deepcopy(mesh.elem)

    # set the z locations
    for k in range(nz):
        for i in range(nel2d):
            iel = i + nel2d*k
            # replace the position arrays with 3D ones (with z empty)
            mesh3d.elem[iel].pos = np.zeros((3, 2, 2, 2))
            mesh3d.elem[iel].pos[:, :1, :, :] = mesh.elem[i].pos
            mesh3d.elem[iel].pos[:, 1:2, :, :] = mesh.elem[i].pos
            
            # fill in the z location
            z1 = zmin + k/nz*(zmax-zmin)
            z2 = zmin + (k+1)/nz*(zmax-zmin)
            mesh3d.elem[iel].pos[2, 0, :, :] = z1
            mesh3d.elem[iel].pos[2, 1, :, :] = z2

    # fix the internal boundary conditions (even though it's probably useless)
    # the end boundary conditions will be overwritten later with the proper ones
    for (iel, el) in enumerate(mesh3d.elem):
        for ibc in range(nbc):
            el.bcs[ibc, 4][0] = 'E'
            el.bcs[ibc, 4][1] = iel+1
            el.bcs[ibc, 4][2] = 5
            el.bcs[ibc, 4][3] = iel+nel2d+1
            el.bcs[ibc, 4][4] = 6
            el.bcs[ibc, 5][0] = 'E'
            el.bcs[ibc, 5][1] = iel+1
            el.bcs[ibc, 5][2] = 6
            el.bcs[ibc, 5][3] = iel-nel2d+1
            el.bcs[ibc, 5][4] = 5
        
    # now fix the end boundary conditions
    # face 5 is at zmin and face 6 is at zmax (with Nek indexing, corresponding to 4 and 5 in Python)
    for i in range(nel2d):
        for ibc in range(nbc):
            i1 = i+(nz-1)*nel2d  # index of the face on the zmax side
            mesh3d.elem[i].bcs[ibc, 4][0] = bc1
            mesh3d.elem[i].bcs[ibc, 4][1] = i+1
            mesh3d.elem[i].bcs[ibc, 4][2] = 5
            mesh3d.elem[i1].bcs[ibc, 5][0] = bc2
            mesh3d.elem[i1].bcs[ibc, 5][1] = i1+1
            mesh3d.elem[i1].bcs[ibc, 5][2] = 6
            # fix the matching faces for the periodic conditions
            if bc1 == 'P':
                mesh3d.elem[i].bcs[ibc, 4][3] = i1+1
                mesh3d.elem[i].bcs[ibc, 4][4] = 6
            if bc2 == 'P':
                mesh3d.elem[i1].bcs[ibc, 5][3] = i+1
                mesh3d.elem[i1].bcs[ibc, 5][4] = 5

    # return the extruded mesh
    return mesh3d
