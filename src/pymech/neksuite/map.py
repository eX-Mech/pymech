import struct

import numpy as np

from ..log import logger


def readma2(fname):
    """A function for reading binary map files (``*.ma2``) for nek5000.

    The map file comtains, for each element in the mesh, the id of the MPI rank that owns it
    followed by the ids of the vertices of the element in the global address space.

    The partitioning is determined by generating the undirected graph formed by the mesh, then
    repeatedly computing and deviding the graph using the Fiedler vector of the graph Laplacian.

    Parameters
    ----------
    fname : str
            file name

    Returns
    -------
    cell  : 2d int ndarray
            list of the vertices for each element of the mesh (in global address space)
    procmap  : 1d int ndarray
            processor map (0-based) indicating ownership for each element of the mesh
    """
    infile = open(fname, "rb")
    #
    # read header
    header = infile.read(132).split()
    nel = int(header[1])

    # number of active elements (nrank - noutflow)
    # nactive = int(header[2])

    # total number of element vertices in the mesh
    # npts  = (2**ldim)*nel
    npts = int(header[5])

    # number of unique element vertices in the mesh
    # nrank = int(header[6])

    # number of points on outflow boundaries ('o  ')
    # noutflow = int(header[7])

    # NOTE: these values can be computed from the others but are included in the
    # header
    # number of levels in the binary partition tree
    # depth = log2(nel)
    # depth = int(header[3])

    # maximum number of elements in partition tree of depth d
    # d2    = 2**d
    # d2 = int(header[4])
    # always double precision
    wdsz = 4
    inttype = "i"

    # detect endianness
    etagb = infile.read(4)
    etagL = struct.unpack("<f", etagb)[0]
    etagL = int(etagL * 1e5) / 1e5
    etagB = struct.unpack(">f", etagb)[0]
    etagB = int(etagB * 1e5) / 1e5
    if etagL == 6.54321:
        logger.debug("Reading little-endian file\n")
        emode = "<"
        # endian = "little"
    elif etagB == 6.54321:
        logger.debug("Reading big-endian file\n")
        emode = ">"
        # endian = "big"
    else:
        raise ValueError("Could not interpret endianness")

    # read the entire contents of the file
    # for each element, there are nvert vertices and a processor id
    nvert = int(npts / nel)  # 2**ldim
    buf = infile.read((nvert + 1) * wdsz * nel)

    # processor map (0-based)
    procmap = np.empty((nel,))
    # list of vertices for each element (in global address space)
    cell = np.empty((nel, nvert))
    for iel in range(nel):
        fi = np.frombuffer(
            buf,
            dtype=emode + inttype,
            count=nvert + 1,
            offset=(nvert + 1) * wdsz * iel,
        )
        procmap[iel] = fi[0]
        cell[iel, :] = fi[1:]
    # close file
    infile.close()
    #
    # output
    return cell, procmap
