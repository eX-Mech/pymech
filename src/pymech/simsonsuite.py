"""Module for reading and writing SIMSON files"""

import struct

import numpy as np

from pymech.core import HexaData

__all__ = ("readdns", "readplane")


# ==============================================================================
def readdns(fname):
    """A function for reading binary data from the SIMSON binary format

    Parameters
    ----------
    fname : str
            file name
    """
    #
    try:
        infile = open(fname, "rb")
    except OSError as e:
        print(f"I/O error ({e.errno}): {e.strerror}")
        return -1
    #
    # ---------------------------------------------------------------------------
    # READ HEADER
    # ---------------------------------------------------------------------------
    wdsz = 8
    realtype = "d"
    #
    # identify endian encoding (and number passive scalar)
    etagb = infile.read(4)
    etagL = struct.unpack("<i", etagb)[0]
    etagB = struct.unpack(">i", etagb)[0]
    # 1644 = 44 + (2*8) * 100 means a maximum of 100 passive scalars
    if (etagL >= 44) & (etagL <= 1644):
        # print('Reading little-endian file\n')
        emode = "<"
        nscal = int((etagL - 44) / (2 * wdsz))
    elif (etagB >= 44) & (etagB <= 1644):
        # print('Reading big-endian file\n')
        emode = ">"
        nscal = int((etagB - 44) / (2 * wdsz))
    else:
        print("ERROR: could not initerpret endianness")
        return -3
    #
    # Reynolds Number
    Re = infile.read(wdsz)
    Re = struct.unpack(emode + realtype, Re)[0]
    #
    # Poiseuille/Couette flag [deprecated]
    PouCou = infile.read(4)
    PouCou = struct.unpack(emode + "i", PouCou)[0]
    #
    # physical box size (here x and z only)
    boxsz = [0 for i in range(3)]
    dum = infile.read(2 * wdsz)
    dum = list(struct.unpack(emode + 2 * realtype, dum))
    boxsz[0] = dum[0]
    boxsz[2] = dum[1]
    #
    # time
    time = infile.read(wdsz)
    time = struct.unpack(emode + realtype, time)[0]
    #
    # dummy variable
    dum = infile.read(wdsz)
    #
    # passive scalar(s) (not tested)
    pr = np.zeros(nscal)
    m = np.zeros(nscal)
    for i in range(nscal):
        pr[i] = infile.read(wdsz)
        pr[i] = struct.unpack(emode + realtype, pr[i])[0]
        m[i] = infile.read(wdsz)
        m[i] = struct.unpack(emode + realtype, m[i])[0]
    #
    # end-of-line
    eol = infile.read(8)
    #
    # box size
    lr1 = infile.read(3 * 4)
    lr1 = list(struct.unpack(emode + 3 * "i", lr1))
    #
    # nfzsym (z-symmetry flag)
    nfzsym = infile.read(4)
    nfzsym = list(struct.unpack(emode + "i", nfzsym))[0]
    #
    # end-of-line
    eol = infile.read(8)
    #
    # compute total number of points per element
    # npel = lr1[0] * lr1[1] * lr1[2]
    #
    # get number of pysical dimensions
    ndim = 2 + (lr1[2] > 1)
    #
    # flow type
    fltype = infile.read(4)
    fltype = struct.unpack(emode + "i", fltype)[0]
    #
    # delta star (and boxsz along y-direction)
    dstar = infile.read(wdsz)
    dstar = struct.unpack(emode + realtype, dstar)[0]
    boxsz[1] = 2 / dstar
    #
    # end-of-line
    eol = infile.read(8)
    #
    # flow-type dependent quantities
    #
    if fltype == -1:
        rlam = infile.read(wdsz)
        rlam = struct.unpack(emode + realtype, rlam)[0]
        eol = infile.read(8)
    if fltype == -2:
        rlam = infile.read(wdsz)
        rlam = struct.unpack(emode + realtype, rlam)[0]
        spanv = infile.read(wdsz)
        spanv = struct.unpack(emode + realtype, spanv)[0]
        eol = infile.read(8)
    if (fltype == 4) or (fltype == 5):
        bstart = infile.read(wdsz)
        bstart = struct.unpack(emode + realtype, bstart)[0]
        blength = infile.read(wdsz)
        blength = struct.unpack(emode + realtype, blength)[0]
        eol = infile.read(8)
    if (fltype >= 4) and (fltype <= 9):
        bstart = infile.read(wdsz)
        bstart = struct.unpack(emode + realtype, bstart)[0]
        blength = infile.read(wdsz)
        blength = struct.unpack(emode + realtype, blength)[0]
        rlam = infile.read(wdsz)
        rlam = struct.unpack(emode + realtype, rlam)[0]
        spanv = infile.read(wdsz)
        spanv = struct.unpack(emode + realtype, spanv)[0]
        eol = infile.read(8)
    if abs(fltype) == 20:
        gr = infile.read(nscal * wdsz)
        gr = struct.unpack(emode + nscal * realtype, gr)
        eol = infile.read(8)
    #
    # get variables
    var = [0 for i in range(5)]
    var[0] = ndim  # position
    var[1] = ndim  # velocity
    var[2] = 0  # pressure is not saved (SIMSON)
    var[3] = 0  # temperature is treated like a scalar (SIMSON)
    var[4] = nscal  # scalars
    #
    # ---------------------------------------------------------------------------
    # READ DATA
    # ---------------------------------------------------------------------------
    #
    # number of points
    #  npel = lr1[0]*lr1[1]*lr1[2]
    #
    # number of points per plane
    nppl = lr1[0] * lr1[2]
    #
    # reading buffer in fourier space
    fou = np.zeros((lr1[2], lr1[1], lr1[0] // 2 + 1)) + 1j * np.zeros(
        (lr1[2], lr1[1], lr1[0] // 2 + 1)
    )
    #
    # initialize data structure
    data = HexaData(ndim, 1, lr1, var)
    data.time = time
    data.wdsz = wdsz
    if emode == "<":
        data.endian = "little"
    elif emode == ">":
        data.endian = "big"
    #
    # generate geometry
    # - x-direction
    dx = boxsz[0] / lr1[0]
    for ix in range(lr1[0]):
        data.elem[0].pos[0, :, :, ix] = dx * ix
    # - y-direction
    dy = np.arccos(-1.0) / (lr1[1] - 1)
    for iy in range(lr1[1]):
        data.elem[0].pos[1, :, iy, :] = boxsz[1] * (1 - np.cos(dy * iy)) / 2
    # - z-direction
    dz = boxsz[2] / lr1[2]
    for iz in range(lr1[2]):
        data.elem[0].pos[2, iz, :, :] = dz * (iz - lr1[2] / 2)
    #
    # read velocity and transform in physical space
    for idim in range(3):
        for iz in range(lr1[2]):
            if iz <= lr1[2] / 2:
                izf = iz
            else:
                izf = lr1[2] // 2 * 3 - (iz + 1)
            for iy in range(lr1[1]):
                fi = infile.read(lr1[0] * wdsz)
                fi = list(struct.unpack(emode + lr1[0] * realtype, fi))
                ip = 0
                for ix in range(int(lr1[0] // 2)):
                    fou[izf, iy, ix] = (fi[ip] + 1j * fi[ip + 1]) * nppl * (-1) ** idim
                    ip += 2
                # end-of-line
                eol = infile.read(8)
        #
        # back to physical space
        data.elem[0].vel[idim, :, :, :] = np.fft.irfft2(fou, (lr1[0], lr1[2]), (2, 0))
    #
    # read scalars and transform in physical space
    for ivar in range(var[4]):
        for iz in range(lr1[2]):
            if iz <= lr1[2] / 2:
                izf = iz
            else:
                izf = lr1[2] // 2 * 3 - (iz + 1)
            for iy in range(lr1[1]):
                fi = infile.read(lr1[0] * wdsz)
                fi = list(struct.unpack(emode + lr1[0] * realtype, fi))
                ip = 0
                for ix in range(int(lr1[0] // 2)):
                    fou[izf, iy, ix] = (fi[ip] + 1j * fi[ip + 1]) * nppl
                    ip += 2
                # end-of-line
                eol = infile.read(8)  # noqa: F841  # required for reading
        #
        # back to physical space
        data.elem[0].scal[ivar, :, :, :] = np.fft.irfft2(fou, (lr1[0], lr1[2]), (2, 0))
    #
    # ---------------------------------------------------------------------------
    # CLOSE FILE
    # ---------------------------------------------------------------------------
    #
    # close file
    infile.close()
    #
    # output
    return data


# ==============================================================================
def readplane(fname):
    """A function for reading binary data from SIMSON's pxyst_ plane files

    .. _pxyst: https://github.com/KTH-Nek5000/SIMSON/tree/master/pxyst

    Parameters
    ----------
    fname : str
            file name
    """
    #
    try:
        infile = open(fname, "rb")
    except OSError as e:
        print(f"I/O error ({e.errno}): {e.strerror}")
        return -1
    #
    # ---------------------------------------------------------------------------
    # READ HEADER
    # ---------------------------------------------------------------------------
    wdsz = 8
    realtype = "d"
    #
    # identify endian encoding (and number passive scalar)
    etagb = infile.read(4)
    etagL = struct.unpack("<i", etagb)[0]
    etagB = struct.unpack(">i", etagb)[0]
    # etagL=1000 for some reason (see rrr.m)
    if etagL <= 1000:
        # print('Reading little-endian file\n')
        endian = "little"
        emode = "<"
        ndim = int(etagL / 4)
    else:
        # print('Reading big-endian file\n')
        endian = "big"
        emode = ">"
        ndim = int(etagB / 4)
    #
    nt = 1
    # nn = np.zeros(ndim)
    nn = []
    for i in range(ndim):
        nnn = infile.read(4)
        nn.append(int(struct.unpack(emode + "i", nnn)[0]))
        nt = nt * nn[i]
    #
    # end-of-line
    eol = infile.read(4)
    #
    if ndim == 1:
        print("reading file %s (%s endian): %d" % (fname, endian, nn[0]))
        x = np.zeros(nn[0])
    elif ndim == 2:
        print("reading file %s (%s endian): %d x %d" % (fname, endian, nn[0], nn[1]))
        x = np.zeros((2, nn[0], nn[1]))
    elif ndim == 3:
        print(
            "reading file %s (%s endian): %d x %d x %d"
            % (fname, endian, nn[0], nn[1], nn[2])
        )
        x = np.zeros((3, nn[0], nn[1], nn[2]))
    else:
        print("ERROR: more than three dimensions")
        return -3
    #
    # ---------------------------------------------------------------------------
    # READ COORDINATES
    # ---------------------------------------------------------------------------
    if ndim == 1:
        # TODO untested
        print("WARNING: reading 1D files was not tested")
        eol = infile.read(4)
        nt1 = struct.unpack(emode + "i", eol)[0] / 8
        dum = infile.read(nt1 * wdsz)
        x = struct.unpack(emode + nt1 * realtype, dum)
        #
        # end-of-line
        eol = infile.read(4)
    elif ndim == 2:
        for i in range(ndim):
            eol = infile.read(4)
            dum = infile.read(nt * wdsz)
            xx = struct.unpack(emode + nt * realtype, dum)
            x[i, :, :] = np.reshape(xx, (nn[0], nn[1]), "F")
            #
            # end-of-line
            eol = infile.read(4)
    elif ndim == 3:
        # TODO untested
        print("WARNING: reading 3D files was not tested")
        for i in range(ndim):
            eol = infile.read(4)
            nt3 = struct.unpack(emode + "i", eol)[0] / 8
            dum = infile.read(nt3 * wdsz)
            xx = struct.unpack(emode + nt3 * realtype, dum)
            x[i, :, :, :] = np.reshape(xx, (nn[0], nn[1], nn[2]), "F")
            #
            # end-of-line
            eol = infile.read(4)
    #
    # ---------------------------------------------------------------------------
    # READ DATA
    # ---------------------------------------------------------------------------
    eol = infile.read(4)
    dum = infile.read(nt * wdsz)
    dd = struct.unpack(emode + nt * realtype, dum)
    if ndim == 1:
        # TODO untested
        d = np.reshape(dd, (nn[0]), "F")
    elif ndim == 2:
        d = np.reshape(dd, (nn[0], nn[1]), "F")
    elif ndim == 3:
        # TODO untested
        d = np.reshape(dd, (nn[0], nn[1], nn[2]), "F")
    #
    # ---------------------------------------------------------------------------
    # CLOSE FILE
    # ---------------------------------------------------------------------------
    #
    # close file
    infile.close()
    #
    # output
    return x, d, nn, ndim
