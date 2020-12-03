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
    # The curved sides will also be extruded, one on each side of each element along nz
    mesh3d.ncurv = 2*nz*mesh.ncurv

    # add extra copies of all elements
    for k in range(nz-1):
        mesh3d.elem = mesh3d.elem + copy.deepcopy(mesh.elem)

    # set the z locations and curvature
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

            # extend curvature and correct it if necessary
            for icurv in range(4):
                curv_type = mesh3d.elem[iel].ccurv[icurv]
                # a 2D element has 4 edges that can be curved, numbered 0-3;
                # the extruded 3D element can have four more (on the other side), numbered 4-7
                mesh3d.elem[iel].ccurv[icurv+4] = curv_type
                mesh3d.elem[iel].curv[icurv+4] = mesh3d.elem[iel].curv[icurv]    # curvature params
                if curv_type == 'm':
                    # in this case the midpoint is given. (x, y) is correct but z should be set to the proper value.
                    mesh3d.elem[iel].curv[icurv][2] = z1
                    mesh3d.elem[iel].curv[icurv+4][2] = z2

    # fix the internal boundary conditions (even though it's probably useless)
    # the end boundary conditions will be overwritten later with the proper ones
    for (iel, el) in enumerate(mesh3d.elem):
        for ibc in range(nbc):
            el.bcs[ibc, 4][0] = 'E'
            el.bcs[ibc, 4][1] = iel+1
            el.bcs[ibc, 4][2] = 5
            el.bcs[ibc, 4][3] = iel-nel2d+1
            el.bcs[ibc, 4][4] = 6
            el.bcs[ibc, 5][0] = 'E'
            el.bcs[ibc, 5][1] = iel+1
            el.bcs[ibc, 5][2] = 6
            el.bcs[ibc, 5][3] = iel+nel2d+1
            el.bcs[ibc, 5][4] = 5
            # update the conditions for side faces
            for iface in range(4):
                el.bcs[ibc, iface][1] = iel+1
                if el.bcs[ibc, iface][0] == 'E':
                    # el.bcs[ibc, 0][1] ought to contain iel+1 once the mesh is valid
                    # but for now it should be off by a factor of nel2d because it is a copy of an element in the first slice
                    offset = iel-el.bcs[ibc, iface][1]+1
                    el.bcs[ibc, iface][3] = el.bcs[ibc, iface][3]+offset
        
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


#==============================================================================
def extrude_split(mesh, z, n, bc1, bc2, Rlim, coord='', ble=0.0, bte=0.0):
    """Extrudes a 2D mesh into a 3D one
    
    Parameters
    ----------
    mesh : exadata
           2D mesh structure to extrude
    z : float
        z values of the boundaries of the intervals of the extruded mesh
    n : int
         number of elements per interval of the extruded mesh
    bc : str
         the boundary condition to use at the ends ('P' is not possible with nonorthogonal extrusion)
    coord : float
            array of length four with the structure (xle,yle,xte,yte) where xle are the x and y coordinates of the leading edge. xte and yte are the counterparts for the trailing edge
    ble : float
          sweep angle (degrees) at the leading edge. A positive ble bends the mesh inward
    bte : float
          sweep angle (degrees) at the trailing edge. A positive bte bends the mesh inward

    Special cases
    -------------
    1) n = [n{1},...,n{i},...,n{N-1}] and z = [z{1},...,z{i},...,z{N}] : The code extrudes the mesh between z{1} and z{N} with n{i} elements in the interval defined by z{i} and z{i+1} (len(n)=len(z)-1)
    1) n = [''] and z = [zmin,zmax] : The code extrudes the mesh between zmin and zmax with the normalized (between 0 and 1) point distribution from z.txt
    2) n = [dz0,s] and z = [zmin,zmax] : The code extrudes the mesh between zmin and zmax with a geometric point distribution defined by the initial spacing dz0 and inflation ratio s
    3) coord = [''] : xle/xte are set to the minimum/maximum x values of the mesh and yle/yte are set to the minimum/maximum y values of the mesh
    4) ble = 0.0 and bte = 0.0 : orthogonal extrusion
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

    if all(v is '' for v in n) and not any(v is '' for v in z) and len(z)==2:
        print('Reading z-point distribution (scaled from 0 to 1) from z.txt file')
        f = open('z.txt')
        lines = f.readlines()
        f.close
        nz = len(lines)-1
        z1 = np.zeros((nz,1))
        z2 = np.zeros((nz,1))
        for k in range(0,nz):
            z1[k] = z[0]+float(lines[k])*(z[1]-z[0])
        for k in range(1,nz+1):
            z2[k-1] = z[0]+float(lines[k])*(z[1]-z[0])
        if len(z1) < 1:
            logger.critical('At least two points are necessary in the z.txt file')
            return -5
        elif len(z1) == 1 or np.amax(np.diff(np.diff(np.append(z1,z2[nz-1])))) < 1e-15:
            print('Uniform extrusion')
        else:
            print('Nonuniform extrusion') 
        mode = 0
    elif not any(v is '' for v in n) and not any(v is '' for v in z) and len(n)==len(z)-1:
        print('Extrusion based on number of elements per interval and z coordinates of intervals')
        nz = np.sum(n)
        z1 = np.zeros((nz,1))
        z2 = np.zeros((nz,1))
        if len(z) < 2:
           logger.critical('z should contain at least two points')
           return -6
        elif len(z) == 2 or np.amax(np.diff(np.divide(np.diff(z),n))) < 1e-15:
            print('Uniform extrusion')
        else:
            print('Nonuniform extrusion')
        it = 0
        for kk in range(len(n)):
            for k in range(0,n[kk]):
                z1[it] = z[kk] + k/float(n[kk])*(z[kk+1]-z[kk])
                z2[it] = z[kk] + (k+1)/float(n[kk])*(z[kk+1]-z[kk])
                it += 1
        mode = 1
    elif not any(v is '' for v in n) and not any(v is '' for v in z) and len(n)==2 and len(z)==2:
        print('Extrusion based on initial spacing and inflation ratio')
        print('Nonuniform extrusion')
        dz0 = n[0]
        s = n[1]
        if s == 1:
            logger.critical('Inflation ratio cannot be 1')
            return -7
        nz = round(np.log(1.0-(z[1]-z[0])/dz0*(1.0-s))/np.log(s))
        z1 = np.zeros((nz,1))
        z2 = np.zeros((nz,1))
        z1[0] = z[0]
        z2[0] = z[0]+dz0
        for k in range(1,nz):
            z1[k] = z1[k-1]+dz0*s**(k-1)
            z2[k] = z2[k-1]+dz0*s**k
        if z2[-1]!=z[1]:
            print('Extrusion exceeds or does not reach requested zmax. Reducing or increasing z2 to zmax')
            z2[-1]=z[1]
            if abs(z1[-1])>=abs(z[1]):
                print('z1 needs to be reduced to be smaller than zmax')
                z1[-1]=0.5*(z1[-2]+z2[-1])
        mode = 2
    else: 
        logger.critical('Inconsistent n and z values')
        return -8
    if nz % 4 != 0:
        logger.critical('Inconsistent elements to extrude: nz must be divided by 4')
 
    # copy the structure and make it 3D
    mesh3d = copy.deepcopy(mesh)
    mesh3d.lr1 = [2, 2, 2]
    mesh3d.var = [3, 0, 0, 0, 0]  # remove anything that isn't geometry for now
    nel2d = mesh.nel
    nel3d = nel2d*6*int(nz/4)
    nbc = mesh.nbc
    mesh3d.nel = nel3d
    mesh3d.ndim = 3
    # The curved sides will also be extruded, one on each side of each element along nz
    mesh3d.ncurv = 2*nz*mesh.ncurv

    # add extra copies of all elements
    for k in range(6*int(nz/4)-1):
        mesh3d.elem = mesh3d.elem + copy.deepcopy(mesh.elem)

    # Change the x and y coordinates if nonorthogonal extrusion
    x1 = np.zeros((nz,2,2))
    x2 = np.zeros((nz,2,2))
    y1 = np.zeros((nz,2,2))
    y2 = np.zeros((nz,2,2))

    if all(v is '' for v in coord):
        print('Leading/trailing edge is set to the minimum/maximum x,y-coordinates in the mesh')
        xvec = np.zeros((4*nel2d,1))
        yvec = np.zeros((4*nel2d,1))
        it = 0
        for iel in range(nel2d):
                for i in range(2):
                    for j in range(2):
                        xvec[it] = mesh.elem[iel].pos[0,0,j,i]
                        yvec[it] = mesh.elem[iel].pos[1,0,j,i]
                        it += 1
        indxmin = np.argmin(xvec)
        indxmax = np.argmax(xvec)
        indymin = np.argmin(yvec)
        indymax = np.argmax(yvec)
        xle = xvec[indxmin]
        yle = yvec[indymin]
        xte = xvec[indxmax]
        yte = yvec[indymax]
        zle = z[0]
        zte = z[0]
    elif not any(v is '' for v in coord):
        print('Leading/trailing edge is set to the x,y-coordinates from the user`s definition')
        xle = coord[0]
        yle = coord[1]
        xte = coord[2]
        yte = coord[3]
        zle = z[0]
        zte = z[0]
    else:
        logger.critical('The x and y positions of the leading and trailing edges are inconsistent')
        return -9

    ble = np.deg2rad(ble)
    bte = np.deg2rad(bte)

    if ble == 0 and bte == 0:
        print('Orthogonal extrusion')
        ort = 1
    else:
        print('Nonorthogonal extrusion')
        ort = 0
        xa = xle+np.tan(ble)*(xte-xle+np.tan(bte)*zte-np.tan(bte)*zle)/(np.tan(ble)+np.tan(bte))
        ya = 0.5*(yle+yte)
        za = (xte-xle+np.tan(bte)*zte+np.tan(ble)*zle)/(np.tan(ble)+np.tan(bte))  

    # set the z locations and curvature
    for k in range(0, nz, 4):
        for i in range(nel2d):
            iel = 6*(i + nel2d*int(k/4))
            for iell in range(6):
                mesh3d.elem[iel+iell] = copy.deepcopy(mesh.elem[i])
                mesh3d.elem[iel+iell].pos = np.zeros((3, 2, 2, 2))
            
            xvec = np.zeros((2,2))
            yvec = np.zeros((2,2))
            rvec = np.zeros((2,2))
            index_low = np.zeros((2,2), dtype=int)
            index_high = np.zeros((2,2), dtype=int)
            iindex_low=0
            iindex_high=0
            for ii in range(2):
                for jj in range(2):
                    xvec[jj,ii] = mesh.elem[i].pos[0,0,jj,ii]
                    yvec[jj,ii] = mesh.elem[i].pos[1,0,jj,ii]
                    rvec[jj,ii] = ((yvec[jj,ii]**2)**0.5)/Rlim - 1.0
                    if rvec[jj,ii] <= 0.0:
                        index_low[iindex_low,:]=[jj,ii]
                        iindex_low += 1
                    else:
                        index_high[iindex_high,:]=[jj,ii]
                        iindex_high += 1
            if (iindex_low != 2) or (iindex_high != 2):
                logger.critical('Mid element not consistent. Criteria must divide elements with 2 points on each side.')
                return -11
            
            poslow=mesh.elem[i].pos[:,0, index_low[:,0], index_low[:,1]]
            posmid=0.5*(mesh.elem[i].pos[:,0, index_low[:,0], index_low[:,1]] + mesh.elem[i].pos[:,0, index_high[:,0], index_high[:,1]])
            poshigh=mesh.elem[i].pos[:,0, index_high[:,0], index_high[:,1]]
            
            mesh3d.elem[iel].pos[:, 0, index_low[:,0], index_low[:,1]] = poslow
            mesh3d.elem[iel].pos[:, 1, index_low[:,0], index_low[:,1]] = poslow
            mesh3d.elem[iel].pos[:, 0, index_high[:,0], index_high[:,1]] = posmid
            mesh3d.elem[iel].pos[:, 1, index_high[:,0], index_high[:,1]] = posmid
            
            mesh3d.elem[iel+1].pos[:, 0, index_low[:,0], index_low[:,1]] = posmid
            mesh3d.elem[iel+1].pos[:, 1, index_low[:,0], index_low[:,1]] = posmid
            mesh3d.elem[iel+1].pos[:, 0, index_high[:,0], index_high[:,1]] = poshigh
            mesh3d.elem[iel+1].pos[:, 1, index_high[:,0], index_high[:,1]] = poshigh
            
            mesh3d.elem[iel+2].pos[:, 0, index_low[:,0], index_low[:,1]] = poslow
            mesh3d.elem[iel+2].pos[:, 1, index_low[:,0], index_low[:,1]] = poslow
            mesh3d.elem[iel+2].pos[:, 0, index_high[:,0], index_high[:,1]] = posmid
            mesh3d.elem[iel+2].pos[:, 1, index_high[:,0], index_high[:,1]] = poshigh
            
            mesh3d.elem[iel+3].pos[:, 0, index_low[:,0], index_low[:,1]] = poslow
            mesh3d.elem[iel+3].pos[:, 1, index_low[:,0], index_low[:,1]] = poslow
            mesh3d.elem[iel+3].pos[:, 0, index_high[:,0], index_high[:,1]] = poshigh
            mesh3d.elem[iel+3].pos[:, 1, index_high[:,0], index_high[:,1]] = posmid
            
            mesh3d.elem[iel+4].pos[:, 0, index_low[:,0], index_low[:,1]] = posmid
            mesh3d.elem[iel+4].pos[:, 1, index_low[:,0], index_low[:,1]] = posmid
            mesh3d.elem[iel+4].pos[:, 0, index_high[:,0], index_high[:,1]] = poshigh
            mesh3d.elem[iel+4].pos[:, 1, index_high[:,0], index_high[:,1]] = poshigh
            
            mesh3d.elem[iel+5].pos[:, 0, index_low[:,0], index_low[:,1]] = poslow
            mesh3d.elem[iel+5].pos[:, 1, index_low[:,0], index_low[:,1]] = poslow
            mesh3d.elem[iel+5].pos[:, 0, index_high[:,0], index_high[:,1]] = posmid
            mesh3d.elem[iel+5].pos[:, 1, index_high[:,0], index_high[:,1]] = posmid
            
            # replace the position arrays with 3D ones (with z empty)
#            mesh3d.elem[iel].pos = np.zeros((3, 2, 2, 2))
#            mesh3d.elem[iel].pos[:, :1, :, :] = mesh.elem[i].pos
#            mesh3d.elem[iel].pos[:, 1:2, :, :] = mesh.elem[i].pos

            if ort == 0:
                # scale and translate x and y coordinates if nonorthogonal extrusion
                # Face 1
                x1[k] = xa+(xa-mesh3d.elem[iel].pos[0, :1, :, :])/(za-z1[0])*(z1[k]-za)
                y1[k] = ya+(ya-mesh3d.elem[iel].pos[1, :1, :, :])/(za-z1[0])*(z1[k]-za)

                mesh3d.elem[iel].pos[0, :1, :, :] = x1[k]
                mesh3d.elem[iel].pos[1, :1, :, :] = y1[k]

                # Face 2
                x2[k] = xa+(xa-mesh3d.elem[iel].pos[0, 1:2, :, :])/(za-z1[0])*(z2[k]-za)
                y2[k] = ya+(ya-mesh3d.elem[iel].pos[1, 1:2, :, :])/(za-z1[0])*(z2[k]-za)

                mesh3d.elem[iel].pos[0, 1:2, :, :] = x2[k]
                mesh3d.elem[iel].pos[1, 1:2, :, :] = y2[k]

            # fill in the z location
            mesh3d.elem[iel].pos[2, 0, :, :] = z1[k]
            mesh3d.elem[iel].pos[2, 1, :, :] = z2[k]
            
            mesh3d.elem[iel+1].pos[2, 0, :, :] = z1[k]
            mesh3d.elem[iel+1].pos[2, 1, index_low[0,0], index_low[0,1]] = z2[k]
            mesh3d.elem[iel+1].pos[2, 1, index_low[1,0], index_low[1,1]] = z2[k]
            mesh3d.elem[iel+1].pos[2, 1, index_high[0,0], index_high[0,1]] = z2[k+1]
            mesh3d.elem[iel+1].pos[2, 1, index_high[1,0], index_high[1,1]] = z2[k+1]
            
            mesh3d.elem[iel+2].pos[2, 0, :, :] = z1[k+1]
            mesh3d.elem[iel+2].pos[2, 1, :, :] = z2[k+1]
            
            mesh3d.elem[iel+3].pos[2, 0, :, :] = z1[k+2]
            mesh3d.elem[iel+3].pos[2, 1, :, :] = z2[k+2]
            
            mesh3d.elem[iel+4].pos[2, 1, :, :] = z2[k+3]
            mesh3d.elem[iel+4].pos[2, 0, index_low[0,0], index_low[0,1]] = z1[k+3]
            mesh3d.elem[iel+4].pos[2, 0, index_low[1,0], index_low[1,1]] = z1[k+3]
            mesh3d.elem[iel+4].pos[2, 0, index_high[0,0], index_high[0,1]] = z1[k+2]
            mesh3d.elem[iel+4].pos[2, 0, index_high[1,0], index_high[1,1]] = z1[k+2]
            
            mesh3d.elem[iel+5].pos[2, 0, :, :] = z1[k+3]
            mesh3d.elem[iel+5].pos[2, 1, :, :] = z2[k+3]

            # extend curvature and correct it if necessary (works only for curv_type == 'C', if all elements have curvature)
            if (index_low[0,:]==[0,0]).all():
                if (index_high[0,:]==[0,1]).all():
                    mesh3d.elem[iel].curv[1]=0.5*(mesh.elem[i].curv[1]-mesh.elem[i].curv[3])
                    mesh3d.elem[iel+1].curv[3]=0.5*(mesh.elem[i].curv[3]-mesh.elem[i].curv[1])
                else:
                    mesh3d.elem[iel].curv[2]=0.5*(mesh.elem[i].curv[2]-mesh.elem[i].curv[0])
                    mesh3d.elem[iel+1].curv[0]=0.5*(mesh.elem[i].curv[0]-mesh.elem[i].curv[2])
            elif (index_low[0,:]==[1,0]).all():
                mesh3d.elem[iel].curv[0]=0.5*(mesh.elem[i].curv[0]-mesh.elem[i].curv[2])
                mesh3d.elem[iel+1].curv[2]=0.5*(mesh.elem[i].curv[2]-mesh.elem[i].curv[0])
            elif (index_low[0,:]==[0,1]).all():
                mesh3d.elem[iel].curv[3]=0.5*(mesh.elem[i].curv[3]-mesh.elem[i].curv[0])
                mesh3d.elem[iel+1].curv[1]=0.5*(mesh.elem[i].curv[1]-mesh.elem[i].curv[3])
            for icurv in range(4):
                # a 2D element has 4 edges that can be curved, numbered 0-3;
                # the extruded 3D element can have four more (on the other side), numbered 4-7
                for iell in range(6):
                    mesh3d.elem[iel+iell].ccurv[icurv+4] = mesh3d.elem[iel].ccurv[icurv]
                    mesh3d.elem[iel+iell].curv[icurv+4] = mesh3d.elem[iel+iell].curv[icurv]    # curvature params
                    
            mesh3d.elem[iel+2].curv[0:4]=mesh3d.elem[iel].curv[0:4]
            mesh3d.elem[iel+3].curv[4:8]=mesh3d.elem[iel].curv[0:4]
            mesh3d.elem[iel+5].curv=mesh3d.elem[iel].curv
            mesh3d.elem[iel+4].curv=mesh3d.elem[iel+1].curv
            mesh3d.elem[iel+2].ccurv[0:4]=mesh3d.elem[iel].ccurv[0:4]
            mesh3d.elem[iel+3].ccurv[4:8]=mesh3d.elem[iel].ccurv[0:4]
            mesh3d.elem[iel+5].ccurv=mesh3d.elem[iel].ccurv
            mesh3d.elem[iel+4].ccurv=mesh3d.elem[iel+1].ccurv
            
            for icurv in range(4):
                curv_type = mesh3d.elem[iel].ccurv[icurv]
                if curv_type == 'm':
                    # x and y are correct if orthogonal extrusion but z should be set to the proper value.
                    mesh3d.elem[iel].curv[icurv][2] = z1[k]
                    mesh3d.elem[iel].curv[icurv+4][2] = z2[k]
                    if ort == 0:
                        # scale and translate x and y coordinates if nonorthogonal extrusion
                        mesh3d.elem[iel].curv[icurv][0]   = xa+(xa-mesh3d.elem[iel].curv[icurv][0])/(za-z1[0])*(z1[k]-za)
                        mesh3d.elem[iel].curv[icurv+4][0] = xa+(xa-mesh3d.elem[iel].curv[icurv+4][0])/(za-z1[0])*(z2[k]-za)
                        mesh3d.elem[iel].curv[icurv][1]   = ya+(ya-mesh3d.elem[iel].curv[icurv][1])/(za-z1[0])*(z1[k]-za)
                        mesh3d.elem[iel].curv[icurv+4][1] = ya+(ya-mesh3d.elem[iel].curv[icurv+4][1])/(za-z1[0])*(z2[k]-za)
                elif curv_type == 'C':
                    if ort == 0:
                        # scaling the radius of the circle
                        mesh3d.elem[iel].curv[icurv][0]   = mesh3d.elem[iel].curv[icurv][0]*abs((z1[k]-za)/(za-z1[0]))
                        mesh3d.elem[iel].curv[icurv+4][0] = mesh3d.elem[iel].curv[icurv+4][0]*abs((z2[k]-za)/(za-z1[0]))
                elif curv_type == 's':
                    if ort == 0:
                        # scaling the radius of the sphere
                        mesh3d.elem[iel].curv[icurv][0]   = mesh3d.elem[iel].curv[icurv][0]*abs((z1[k]-za)/(za-z1[0]))
                        mesh3d.elem[iel].curv[icurv+4][0] = mesh3d.elem[iel].curv[icurv+4][0]*abs((z2[k]-za)/(za-z1[0])) 
                    # x and y are correct if orthogonal extrusion but z should be set to the proper value. 
                    mesh3d.elem[iel].curv[icurv][3] = z1[k]
                    mesh3d.elem[iel].curv[icurv+4][3] = z2[k]
                   
                    if ort == 0:
                        # scale and translate x and y coordinates if nonorthogonal extrusion
                        mesh3d.elem[iel].curv[icurv][1]   = xa+(xa-mesh3d.elem[iel].curv[icurv][1])/(za-z1[0])*(z1[k]-za)
                        mesh3d.elem[iel].curv[icurv+4][1] = xa+(xa-mesh3d.elem[iel].curv[icurv+4][1])/(za-z1[0])*(z2[k]-za)
                        mesh3d.elem[iel].curv[icurv][2]   = ya+(ya-mesh3d.elem[iel].curv[icurv][2])/(za-z1[0])*(z1[k]-za)
                        mesh3d.elem[iel].curv[icurv+4][2] = ya+(ya-mesh3d.elem[iel].curv[icurv+4][2])/(za-z1[0])*(z2[k]-za)

    # fix the internal boundary conditions (even though it's probably useless) !!WRONG FIXME
    # the end boundary conditions will be overwritten later with the proper ones
#    for (iel, el) in enumerate(mesh3d.elem):
            for ibc in range(nbc):
                # set the conditions in faces normal to z
                for iell in range(6):
                    mesh3d.elem[iel+iell].bcs[ibc, 4][0] = 'E'
                    mesh3d.elem[iel+iell].bcs[ibc, 4][1] = iel+iell+1
                    mesh3d.elem[iel+iell].bcs[ibc, 4][2] = 5
                    mesh3d.elem[iel+iell].bcs[ibc, 4][4] = 6
                    mesh3d.elem[iel+iell].bcs[ibc, 5][0] = 'E'
                    mesh3d.elem[iel+iell].bcs[ibc, 5][1] = iel+iell+1
                    mesh3d.elem[iel+iell].bcs[ibc, 5][2] = 6
                    mesh3d.elem[iel+iell].bcs[ibc, 5][4] = 5
                
                mesh3d.elem[iel].bcs[ibc, 4][3] = iel+1-nel2d
                mesh3d.elem[iel].bcs[ibc, 5][3] = iel+1+2
                mesh3d.elem[iel+1].bcs[ibc, 4][3] = iel+1-nel2d
                mesh3d.elem[iel+1].bcs[ibc, 5][3] = iel+1+2
                mesh3d.elem[iel+2].bcs[ibc, 4][3] = iel+1
                mesh3d.elem[iel+2].bcs[ibc, 5][3] = iel+1+3
                mesh3d.elem[iel+3].bcs[ibc, 4][3] = iel+1+2
                mesh3d.elem[iel+3].bcs[ibc, 5][3] = iel+1+5
                mesh3d.elem[iel+4].bcs[ibc, 4][3] = iel+1+3
                mesh3d.elem[iel+4].bcs[ibc, 5][3] = iel+1+nel2d
                mesh3d.elem[iel+5].bcs[ibc, 4][3] = iel+1+3
                mesh3d.elem[iel+5].bcs[ibc, 5][3] = iel+1+nel2d
                # update the conditions for side faces
                for iface in range(4):
                    mesh3d.elem[iel].bcs[ibc, iface][1] = iel+1
                    if mesh3d.elem[iel].bcs[ibc, iface][0] == 'E':
                        # el.bcs[ibc, 0][1] ought to contain iel+1 once the mesh is valid
                        # but for now it should be off by a factor of nel2d because it is a copy of an element in the first slice
                        offset = iel-mesh3d.elem[iel].bcs[ibc, iface][1]+1
                        mesh3d.elem[iel].bcs[ibc, iface][1] = iel+1
                        mesh3d.elem[iel].bcs[ibc, iface][3] = mesh3d.elem[iel].bcs[ibc, iface][3]+offset
             
    # now fix the end boundary conditions 
    #!!WRONG FIXME if original element is part of boundary (creates the same wall or velocity boundary condition on elements i to i+5, while the right would be i,i+2,i+3 and i+5 (internal) or i+1 and i+4 (external))
    # face 5 is at zmin and face 6 is at zmax (with Nek indexing, corresponding to 4 and 5 in Python)
    for i in range(0,6*nel2d,6):
        for ibc in range(nbc):
            i1 = i+nel3d-6*nel2d+5
            mesh3d.elem[i].bcs[ibc, 4][0] = bc1
            mesh3d.elem[i].bcs[ibc, 4][1] = i+1
            mesh3d.elem[i].bcs[ibc, 4][2] = 5
            mesh3d.elem[i+1].bcs[ibc, 4][0] = bc1
            mesh3d.elem[i+1].bcs[ibc, 4][1] = i+1+1
            mesh3d.elem[i+1].bcs[ibc, 4][2] = 5
            mesh3d.elem[i1].bcs[ibc, 5][0] = bc2
            mesh3d.elem[i1].bcs[ibc, 5][1] = i1+1
            mesh3d.elem[i1].bcs[ibc, 5][2] = 6
            mesh3d.elem[i1-1].bcs[ibc, 5][0] = bc2
            mesh3d.elem[i1-1].bcs[ibc, 5][1] = i1-1+1
            mesh3d.elem[i1-1].bcs[ibc, 5][2] = 6
            # fix the matching faces for the periodic conditions
            if bc1 == 'P':
                mesh3d.elem[i].bcs[ibc, 4][3] = i1+1
                mesh3d.elem[i].bcs[ibc, 4][4] = 6
                mesh3d.elem[i+1].bcs[ibc, 4][3] = i1-1+1
                mesh3d.elem[i+1].bcs[ibc, 4][4] = 6
            if bc2 == 'P':
                mesh3d.elem[i1].bcs[ibc, 5][3] = i+1
                mesh3d.elem[i1].bcs[ibc, 5][4] = 5
                mesh3d.elem[i1-1].bcs[ibc, 5][3] = i+1+1
                mesh3d.elem[i1-1].bcs[ibc, 5][4] = 5
    
    # FIND THE CURVED ELEMENTS
    ncurv = 0
    for el in mesh3d.elem:
        for iedge in range(12):
            if el.ccurv[iedge] != '':
                ncurv = ncurv+1
    mesh3d.ncurv = ncurv
    return mesh3d
