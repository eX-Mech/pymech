import struct

import numpy as np
from .log import logger


def readma2(fname,ldim):
    """A function for reading binary map files (*.ma2) for nek5000.

    The map file comtains, for each element in the mesh, the id of the MPI rank that owns it 
    followed by the ids of the vertices of the element in the global address space.

    The partitioning is determined by generating the undirected graph formed by the mesh, then 
    repeatedly computing and deviding the graph using the Fiedler vector of the graph Laplacian. 

    Parameters
    ----------
    fname : str
            file name
    ldim  : int
            dimension of the mesh
    """
    try:
        infile = open(fname, "rb")
    except OSError as e:
        logger.critical(f"I/O error ({e.errno}): {e.strerror}")
        return -1
    #
    # read header
    header   = infile.read(132).split()
    nel      = int(header[1])
    nactive  = int(header[2])
    depth    = int(header[3])
    d2       = int(header[4])
    npts     = int(header[5])
    nrnk     = int(header[6])
    noutflow = int(header[7])
    # always double precision
    wdsz    = 4
    inttype = 'i'

    # detect endianness
    etagb = infile.read(4)
    etagL = struct.unpack("<f", etagb)[0]
    etagL = int(etagL * 1e5) / 1e5
    etagB = struct.unpack(">f", etagb)[0]
    etagB = int(etagB * 1e5) / 1e5
    if etagL == 6.54321:
        logger.debug("Reading little-endian file\n")
        emode = "<"
        endian = "little"
    elif etagB == 6.54321:
        logger.debug("Reading big-endian file\n")
        emode = ">"
        endian = "big"
    else:
        logger.error("Could not interpret endianness")
        return -3

    # read the entire contents of the file
    # for each element, there are nvert vertices and a processor id
    nvert = 2**ldim
    buf = infile.read((nvert + 1) * wdsz * nel)
    
    # processor map (0-based)
    pmap = np.empty((nel,))
    # list of vertices for each element (in global address space)
    cell = np.empty((nel,nvert))
    for iel in range(nel):
        fi = np.frombuffer(
                buf,
                dtype  = emode + inttype,
                count  = nvert + 1,
                offset = (nvert + 1) * wdsz * iel
            )
        pmap[iel] = fi[0]
        cell[iel,:] = fi[1:]
    # close file
    infile.close()
    #
    # output
    return pmap, cell
