from __future__ import annotations

import struct
import sys

import numpy as np

from pymech.core import HexaData
from pymech.log import logger


def readrea(fname):
    """A function for reading .rea files for nek5000

    Parameters
    ----------
    fname : str
            file name
    """
    #
    infile = open(fname)
    #
    # ---------------------------------------------------------------------------
    # count the number of boundary conditions
    # (it's too dangerous to infer it from the header)
    # ---------------------------------------------------------------------------
    #
    nbc = 0
    for line in infile:
        line_split = line.split()
        if "BOUNDARY" in line_split[2:] and "NO" not in line_split:
            nbc = nbc + 1

    infile.seek(0)
    #
    # ---------------------------------------------------------------------------
    # READ HEADER (2 lines) + ndim + number of parameters
    # ---------------------------------------------------------------------------
    #
    infile.readline()
    infile.readline()
    ndim = int(infile.readline().split()[0])
    npar = int(infile.readline().split()[0])
    #
    nface = 2 * ndim
    #
    # ---------------------------------------------------------------------------
    # READ parameters
    # ---------------------------------------------------------------------------
    #
    param = np.zeros((npar, 1))
    for ipar in range(npar):
        param[ipar] = float(infile.readline().split()[0])
    #
    # ---------------------------------------------------------------------------
    # skip passive scalars
    # ---------------------------------------------------------------------------
    #
    npscal_data = int(infile.readline().split()[0])
    for ipscal in range(npscal_data):
        infile.readline()
    #
    # ---------------------------------------------------------------------------
    # skip logical switches
    # ---------------------------------------------------------------------------
    #
    nswitch = int(infile.readline().split()[0])
    for iswitch in range(nswitch):
        infile.readline()
    #
    # ---------------------------------------------------------------------------
    # skip XFAC,YFAC,XZERO,YZERO
    # ---------------------------------------------------------------------------
    #
    infile.readline()
    #
    # ---------------------------------------------------------------------------
    # READ MESH
    # ---------------------------------------------------------------------------
    #
    infile.readline()
    nel = int(infile.readline().split()[0])
    #
    # initialize data structure
    lr1 = [2, 2, ndim - 1]
    var = [ndim, 0, 0, 0, 0]
    #
    data = HexaData(ndim, nel, lr1, var, nbc)
    #
    # read geometry
    for iel in range(nel):
        # skip element number and group
        infile.readline()
        for idim in range(var[0] - 1):  # if ndim == 3 do this twice
            for jdim in range(var[0]):
                fi = infile.readline().split()
                data.elem[iel].pos[jdim, idim, 0, 0] = float(fi[0])
                data.elem[iel].pos[jdim, idim, 0, 1] = float(fi[1])
                data.elem[iel].pos[jdim, idim, 1, 1] = float(fi[2])
                data.elem[iel].pos[jdim, idim, 1, 0] = float(fi[3])
    #
    # ---------------------------------------------------------------------------
    # CURVED SIDE DATA
    # ---------------------------------------------------------------------------
    #
    infile.readline()
    ncurved = int(infile.readline().split()[0])
    data.ncurv = ncurved
    for icurved in range(ncurved):
        line = infile.readline()
        if nel < 1e3:
            iedge = int(line[0:3]) - 1
            iel = int(line[3:6]) - 1
            data.elem[iel].curv[iedge][0] = float(line[6:20])
            data.elem[iel].curv[iedge][1] = float(line[20:34])
            data.elem[iel].curv[iedge][2] = float(line[34:48])
            data.elem[iel].curv[iedge][3] = float(line[48:62])
            data.elem[iel].curv[iedge][4] = float(line[62:76])
            data.elem[iel].ccurv[iedge] = line[76:79].split()[0]
        elif nel < 1e6:
            iedge = int(line[0:2]) - 1
            iel = int(line[2:8]) - 1
            data.elem[iel].curv[iedge][0] = float(line[8:22])
            data.elem[iel].curv[iedge][1] = float(line[22:36])
            data.elem[iel].curv[iedge][2] = float(line[36:50])
            data.elem[iel].curv[iedge][3] = float(line[50:64])
            data.elem[iel].curv[iedge][4] = float(line[64:78])
            data.elem[iel].ccurv[iedge] = line[78:81].split()[0]
        else:
            iedge = int(line[0:2]) - 1
            iel = int(line[2:12]) - 1
            data.elem[iel].curv[iedge][0] = float(line[12:26])
            data.elem[iel].curv[iedge][1] = float(line[26:40])
            data.elem[iel].curv[iedge][2] = float(line[40:54])
            data.elem[iel].curv[iedge][3] = float(line[54:68])
            data.elem[iel].curv[iedge][4] = float(line[68:82])
            data.elem[iel].ccurv[iedge] = line[82:85].split()[0]
    #
    # ---------------------------------------------------------------------------
    # BOUNDARY CONDITIONS
    # ---------------------------------------------------------------------------
    #
    infile.readline()  # ***** BOUNDARY CONDITIONS *****
    for ibc in range(nbc):
        infile.readline()  # ***** FLUID   BOUNDARY CONDITIONS ***** [or similar]
        for iel in range(nel):
            for iface in range(nface):
                line = infile.readline()
                if nel < 1e3:
                    data.elem[iel].bcs[ibc, iface][0] = line[1:3].strip()
                    data.elem[iel].bcs[ibc, iface][1] = int(line[4:7])
                    data.elem[iel].bcs[ibc, iface][2] = int(line[7:10])
                    data.elem[iel].bcs[ibc, iface][3] = float(line[10:24])
                    data.elem[iel].bcs[ibc, iface][4] = float(line[24:38])
                    data.elem[iel].bcs[ibc, iface][5] = float(line[38:52])
                    data.elem[iel].bcs[ibc, iface][6] = float(line[52:66])
                    data.elem[iel].bcs[ibc, iface][7] = float(line[66:80])
                elif nel < 1e6:
                    data.elem[iel].bcs[ibc, iface][0] = line[1:3].strip()
                    data.elem[iel].bcs[ibc, iface][1] = iel + 1
                    data.elem[iel].bcs[ibc, iface][2] = iface + 1
                    data.elem[iel].bcs[ibc, iface][3] = float(line[10:24])
                    data.elem[iel].bcs[ibc, iface][4] = float(line[24:38])
                    data.elem[iel].bcs[ibc, iface][5] = float(line[38:52])
                    data.elem[iel].bcs[ibc, iface][6] = float(line[52:66])
                    data.elem[iel].bcs[ibc, iface][7] = float(line[66:80])
                else:
                    data.elem[iel].bcs[ibc, iface][0] = line[1:3].strip()
                    data.elem[iel].bcs[ibc, iface][1] = int(line[4:15])
                    data.elem[iel].bcs[ibc, iface][2] = int(line[15:16])
                    data.elem[iel].bcs[ibc, iface][3] = float(line[16:34])
                    data.elem[iel].bcs[ibc, iface][4] = float(line[34:52])
                    data.elem[iel].bcs[ibc, iface][5] = float(line[52:70])
                    data.elem[iel].bcs[ibc, iface][6] = float(line[70:88])
                    data.elem[iel].bcs[ibc, iface][7] = float(line[88:106])
                # ignore some invalid internal 'E' conditions.
                # They are typically written this way by re2torea and Nek5000 ignores them.
                if (
                    data.elem[iel].bcs[ibc, iface][0] == "E"
                    and data.elem[iel].bcs[ibc, iface][3] == 0.0
                ):
                    data.elem[iel].bcs[ibc, iface][0] = ""
                    for j in range(1, 8):
                        data.elem[iel].bcs[ibc, iface][j] = 0
    #
    # ---------------------------------------------------------------------------
    # FORGET ABOUT WHAT FOLLOWS
    # ---------------------------------------------------------------------------
    #
    #
    # close file
    infile.close()
    #
    # output
    return data


# ==============================================================================
def writerea(fname, data):
    """A function for writing ascii .rea files for nek5000

    Parameters
    ----------
    fname : str
            file name
    data : :class:`pymech.core.HexaData`
            data structure
    """
    #
    outfile = open(fname, "w")
    #
    # ---------------------------------------------------------------------------
    # READ HEADER (2 lines) + ndim + number of parameters
    # ---------------------------------------------------------------------------
    #
    outfile.write("****** PARAMETERS ******\n")
    outfile.write("   2.6000     NEKTON VERSION\n")
    outfile.write(f"   {data.ndim:1d} DIMENSIONAL RUN\n")
    outfile.write("         118 PARAMETERS FOLLOW\n")
    outfile.write("   1.00000     P001: DENSITY\n")
    outfile.write("  -1000.00     P002: VISCOSITY\n")
    outfile.write("   0.00000     P003: BETAG\n")
    outfile.write("   0.00000     P004: GTHETA\n")
    outfile.write("   0.00000     P005: PGRADX\n")
    outfile.write("   0.00000     P006: \n")
    outfile.write("   1.00000     P007: RHOCP\n")
    outfile.write("   1.00000     P008: CONDUCT\n")
    outfile.write("   0.00000     P009: \n")
    outfile.write("   0.00000     P010: FINTIME\n")
    outfile.write("   103         P011: NSTEPS\n")
    outfile.write("  -1.00000E-03 P012: DT\n")
    outfile.write("   0.00000     P013: IOCOMM\n")
    outfile.write("   0.00000     P014: IOTIME\n")
    outfile.write("   10          P015: IOSTEP\n")
    outfile.write("   0.00000     P016: PSSOLVER: 0=default\n")
    outfile.write("   1.00000     P017: \n")
    outfile.write("   0.00000     P018: GRID <0 --> # cells on screen\n")
    outfile.write("   0.00000     P019: INTYPE\n")
    outfile.write("   10.0000     P020: NORDER\n")
    outfile.write("   1.00000E-09 P021: DIVERGENCE\n")
    outfile.write("   1.00000E-09 P022: HELMHOLTZ\n")
    outfile.write("   0.00000     P023: NPSCAL\n")
    outfile.write("   1.00000E-02 P024: TOLREL\n")
    outfile.write("   1.00000E-02 P025: TOLABS\n")
    outfile.write("   1.00000     P026: COURANT/NTAU\n")
    outfile.write("   3.00000     P027: TORDER\n")
    outfile.write("   0.00000     P028: TORDER: mesh velocity (0: p28=p27)\n")
    outfile.write("   0.00000     P029: = magnetic visc if > 0, = -1/Rm if < 0\n")
    outfile.write("   0.00000     P030: > 0 ==> properties set in uservp()\n")
    outfile.write("   0.00000     P031: NPERT: #perturbation modes\n")
    outfile.write("   0.00000     P032: #BCs in re2 file, if > 0\n")
    outfile.write("   0.00000     P033: \n")
    outfile.write("   0.00000     P034: \n")
    outfile.write("   0.00000     P035: \n")
    outfile.write("   0.00000     P036: XMAGNET\n")
    outfile.write("   0.00000     P037: NGRIDS\n")
    outfile.write("   0.00000     P038: NORDER2\n")
    outfile.write("   0.00000     P039: NORDER3\n")
    outfile.write("   0.00000     P040: \n")
    outfile.write("   0.00000     P041: 1-->multiplicattive SEMG\n")
    outfile.write("   0.00000     P042: 0=gmres/1=pcg\n")
    outfile.write("   0.00000     P043: 0=semg/1=schwarz\n")
    outfile.write("   0.00000     P044: 0=E-based/1=A-based prec.\n")
    outfile.write("   0.00000     P045: Relaxation factor for DTFS\n")
    outfile.write("   0.00000     P046: reserved\n")
    outfile.write("   0.00000     P047: vnu: mesh material prop.\n")
    outfile.write("   0.00000     P048: \n")
    outfile.write("   0.00000     P049: \n")
    outfile.write("   0.00000     P050: \n")
    outfile.write("   0.00000     P051: \n")
    outfile.write("   0.00000     P052: IOHIS\n")
    outfile.write("   0.00000     P053: \n")
    outfile.write("   0.00000     P054: fixed flow rate dir: |p54|=1,2,3=x,y,z\n")
    outfile.write("   0.00000     P055: vol.flow rate (p54>0) or Ubar (p54<0)\n")
    outfile.write("   0.00000     P056: \n")
    outfile.write("   0.00000     P057: \n")
    outfile.write("   0.00000     P058: \n")
    outfile.write("   0.00000     P059: !=0 --> full Jac. eval. for each el.\n")
    outfile.write("   0.00000     P060: !=0 --> init. velocity to small nonzero\n")
    outfile.write("   0.00000     P061: \n")
    outfile.write("   0.00000     P062: >0 --> force byte_swap for output\n")
    outfile.write("   8.00000     P063: =8 --> force 8-byte output\n")
    outfile.write("   0.00000     P064: =1 --> perturbation restart\n")
    outfile.write("   1.00000     P065: #iofiles (eg, 0 or 64); <0 --> sep. dirs\n")
    outfile.write("   6.00000     P066: output : <0=ascii, else binary\n")
    outfile.write("   6.00000     P067: restart: <0=ascii, else binary\n")
    outfile.write("   0.00000     P068: iastep: freq for avg_all (0=iostep)\n")
    outfile.write("   0.00000     P069:       : freq of srf dump\n")
    outfile.write("   0.00000     P070: \n")
    outfile.write("   0.00000     P071: \n")
    outfile.write("   0.00000     P072: \n")
    outfile.write("   0.00000     P073: \n")
    outfile.write("   0.00000     P074: \n")
    outfile.write("   0.00000     P075: \n")
    outfile.write("   0.00000     P076: \n")
    outfile.write("   0.00000     P077: \n")
    outfile.write("   0.00000     P078: \n")
    outfile.write("   0.00000     P079: \n")
    outfile.write("   0.00000     P080: \n")
    outfile.write("   0.00000     P081: \n")
    outfile.write("   0.00000     P082: \n")
    outfile.write("   0.00000     P083: \n")
    outfile.write("   0.00000     P084: != 0 --> sets initial timestep if p12>0\n")
    outfile.write("   0.00000     P085: dt retio of p84 !=0, for timesteps>0\n")
    outfile.write("   0.00000     P086: reserved\n")
    outfile.write("   0.00000     P087: \n")
    outfile.write("   0.00000     P088: \n")
    outfile.write("   0.00000     P089: \n")
    outfile.write("   0.00000     P090: \n")
    outfile.write("   0.00000     P091: \n")
    outfile.write("   0.00000     P092: \n")
    outfile.write("   20.0000     P093: Number of previous pressure solns saved\n")
    outfile.write("   9.00000     P094: start projecting velocity after p94 step\n")
    outfile.write("   9.00000     P095: start projecting pressure after p95 step\n")
    outfile.write("   0.00000     P096: \n")
    outfile.write("   0.00000     P097: \n")
    outfile.write("   0.00000     P098: \n")
    outfile.write("   3.00000     P099: dealiasing: <0--> off /3--> old /4-->new\n")
    outfile.write("   0.00000     P100: \n")
    outfile.write("   0.00000     P101: Number of additional modes to filter\n")
    outfile.write("   0.00000     P102: Dump out divergence at each time step\n")
    outfile.write("   0.01000     P103: weight of stabilizing filter\n")
    outfile.write("   0.00000     P104: \n")
    outfile.write("   0.00000     P105: \n")
    outfile.write("   0.00000     P106: \n")
    outfile.write("   0.00000     P107: !=0 --> add h2 array in hmholtz eqn\n")
    outfile.write("   0.00000     P108: \n")
    outfile.write("   0.00000     P109: \n")
    outfile.write("   0.00000     P110: \n")
    outfile.write("   0.00000     P111: \n")
    outfile.write("   0.00000     P112: \n")
    outfile.write("   0.00000     P113: \n")
    outfile.write("   0.00000     P114: \n")
    outfile.write("   0.00000     P115: \n")
    outfile.write("   0.00000     P116: \n")
    outfile.write("   0.00000     P117: \n")
    outfile.write("   0.00000     P118: \n")
    outfile.write("      4  Lines of passive scalar data follows2 CONDUCT, 2RHOCP\n")
    outfile.write(
        "   1.00000        1.00000        1.00000        1.00000        1.00000\n"
    )
    outfile.write("   1.00000        1.00000        1.00000        1.00000\n")
    outfile.write(
        "   1.00000        1.00000        1.00000        1.00000        1.00000\n"
    )
    outfile.write("   1.00000        1.00000        1.00000        1.00000\n")
    outfile.write("         13   LOGICAL SWITCHES FOLLOW\n")
    outfile.write(" T      IFFLOW\n")
    outfile.write(" F      IFHEAT\n")
    outfile.write(" T      IFTRAN\n")
    outfile.write(
        " T F F F F F F F F F F  IFNAV & IFADVC (convection in P.S. fields)\n"
    )
    outfile.write(
        " F F T T T T T T T T T T  IFTMSH (IF mesh for this field is T mesh)\n"
    )
    outfile.write(" F      IFAXIS\n")
    outfile.write(" F      IFSTRS\n")
    outfile.write(" F      IFSPLIT\n")
    outfile.write(" F      IFMGRID\n")
    outfile.write(" F      IFMODEL\n")
    outfile.write(" F      IFKEPS\n")
    outfile.write(" F      IFMVBD\n")
    outfile.write(" F      IFCHAR\n")
    outfile.write(
        "   2.00000       2.00000      -1.00000      -1.00000     XFAC,YFAC,XZERO,YZERO\n"
    )
    #
    # vertex data
    outfile.write(
        "  ***** MESH DATA *****  6 lines are X,Y,Z;X,Y,Z. Columns corners 1-4;5-8\n"
    )
    outfile.write(f"  {data.nel:10d} {data.ndim:1d} {data.nel:10d} NEL,NDIM,NELV\n")
    for iel in range(data.nel):
        outfile.write("           ELEMENT {:10d} [  1a]    GROUP 1\n".format(iel + 1))
        if data.ndim == 2:
            outfile.write(
                "{:14.6e}{:14.6e}{:14.6e}{:14.6e}\n".format(
                    data.elem[iel].pos[0, 0, 0, 0],
                    data.elem[iel].pos[0, 0, 0, -1],
                    data.elem[iel].pos[0, 0, -1, -1],
                    data.elem[iel].pos[0, 0, -1, 0],
                )
            )
            outfile.write(
                "{:14.6e}{:14.6e}{:14.6e}{:14.6e}\n".format(
                    data.elem[iel].pos[1, 0, 0, 0],
                    data.elem[iel].pos[1, 0, 0, -1],
                    data.elem[iel].pos[1, 0, -1, -1],
                    data.elem[iel].pos[1, 0, -1, 0],
                )
            )
        elif data.ndim == 3:
            outfile.write(
                "{:14.6e}{:14.6e}{:14.6e}{:14.6e}\n".format(
                    data.elem[iel].pos[0, 0, 0, 0],
                    data.elem[iel].pos[0, 0, 0, -1],
                    data.elem[iel].pos[0, 0, -1, -1],
                    data.elem[iel].pos[0, 0, -1, 0],
                )
            )
            outfile.write(
                "{:14.6e}{:14.6e}{:14.6e}{:14.6e}\n".format(
                    data.elem[iel].pos[1, 0, 0, 0],
                    data.elem[iel].pos[1, 0, 0, -1],
                    data.elem[iel].pos[1, 0, -1, -1],
                    data.elem[iel].pos[1, 0, -1, 0],
                )
            )
            outfile.write(
                "{:14.6e}{:14.6e}{:14.6e}{:14.6e}\n".format(
                    data.elem[iel].pos[2, 0, 0, 0],
                    data.elem[iel].pos[2, 0, 0, -1],
                    data.elem[iel].pos[2, 0, -1, -1],
                    data.elem[iel].pos[2, 0, -1, 0],
                )
            )
            #
            outfile.write(
                "{:14.6e}{:14.6e}{:14.6e}{:14.6e}\n".format(
                    data.elem[iel].pos[0, -1, 0, 0],
                    data.elem[iel].pos[0, -1, 0, -1],
                    data.elem[iel].pos[0, -1, -1, -1],
                    data.elem[iel].pos[0, -1, -1, 0],
                )
            )
            outfile.write(
                "{:14.6e}{:14.6e}{:14.6e}{:14.6e}\n".format(
                    data.elem[iel].pos[1, -1, 0, 0],
                    data.elem[iel].pos[1, -1, 0, -1],
                    data.elem[iel].pos[1, -1, -1, -1],
                    data.elem[iel].pos[1, -1, -1, 0],
                )
            )
            outfile.write(
                "{:14.6e}{:14.6e}{:14.6e}{:14.6e}\n".format(
                    data.elem[iel].pos[2, -1, 0, 0],
                    data.elem[iel].pos[2, -1, 0, -1],
                    data.elem[iel].pos[2, -1, -1, -1],
                    data.elem[iel].pos[2, -1, -1, 0],
                )
            )
    #
    # curved side data
    outfile.write("  ***** CURVED SIDE DATA *****\n")
    outfile.write(
        f"  {data.ncurv:10d} Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE\n"
    )
    for iel in range(data.nel):
        if data.nel < 1e3:
            #
            for iedge in range(12):
                if data.elem[iel].ccurv[iedge] != "":
                    outfile.write(
                        "{:3d}{:3d}{:14.6e}{:14.6e}{:14.6e}{:14.6e}{:14.6e}{:>2s}\n".format(
                            iedge + 1,
                            iel + 1,
                            data.elem[iel].curv[iedge][0],
                            data.elem[iel].curv[iedge][1],
                            data.elem[iel].curv[iedge][2],
                            data.elem[iel].curv[iedge][3],
                            data.elem[iel].curv[iedge][4],
                            data.elem[iel].ccurv[iedge],
                        )
                    )
        elif data.nel < 1e6:
            #
            for iedge in range(12):
                if data.elem[iel].ccurv[iedge] != "":
                    outfile.write(
                        "{:2d}{:6d}{:14.6e}{:14.6e}{:14.6e}{:14.6e}{:14.6e}{:>2s}\n".format(
                            iedge + 1,
                            iel + 1,
                            data.elem[iel].curv[iedge][0],
                            data.elem[iel].curv[iedge][1],
                            data.elem[iel].curv[iedge][2],
                            data.elem[iel].curv[iedge][3],
                            data.elem[iel].curv[iedge][4],
                            data.elem[iel].ccurv[iedge],
                        )
                    )
        else:
            #
            for iedge in range(12):
                if data.elem[iel].ccurv[iedge] != "":
                    outfile.write(
                        "{:2d}{:10d}{:14.6e}{:14.6e}{:14.6e}{:14.6e}{:14.6e}{:>2s}\n".format(
                            iedge + 1,
                            iel + 1,
                            data.elem[iel].curv[iedge][0],
                            data.elem[iel].curv[iedge][1],
                            data.elem[iel].curv[iedge][2],
                            data.elem[iel].curv[iedge][3],
                            data.elem[iel].curv[iedge][4],
                            data.elem[iel].ccurv[iedge],
                        )
                    )
    #
    # boundary conditions data
    outfile.write("  ***** BOUNDARY CONDITIONS *****\n")
    outfile.write("  ***** FLUID BOUNDARY CONDITIONS *****\n")
    for ibc in range(data.nbc):
        if ibc == 1:
            outfile.write("  ***** THERMAL BOUNDARY CONDITIONS *****\n")
        elif ibc > 1:
            outfile.write(
                "  ***** PASSIVE SCALAR {:4d} BOUNDARY CONDITIONS *****\n".format(
                    ibc - 1
                )
            )
        for iel in range(data.nel):
            for iface in range(2 * data.ndim):
                # if no boundary condition is specified, write 'E' iel iface [...]
                # this is the behaviour of re2torea.
                if data.elem[iel].bcs[ibc, iface][0] == "":
                    data.elem[iel].bcs[ibc, iface][0] = "E"
                    data.elem[iel].bcs[ibc, iface][1] = iel + 1
                    data.elem[iel].bcs[ibc, iface][2] = iface + 1
                if data.nel < 1e3:
                    outfile.write(
                        " {:2s} {:3d}{:3d}{:14.6e}{:14.6e}{:14.6e}{:14.6e}{:14.6e}\n".format(
                            data.elem[iel].bcs[ibc, iface][0],
                            data.elem[iel].bcs[ibc, iface][1],
                            data.elem[iel].bcs[ibc, iface][2],
                            data.elem[iel].bcs[ibc, iface][3],
                            data.elem[iel].bcs[ibc, iface][4],
                            data.elem[iel].bcs[ibc, iface][5],
                            data.elem[iel].bcs[ibc, iface][6],
                            data.elem[iel].bcs[ibc, iface][7],
                        )
                    )
                elif data.nel < 1e6:
                    outfile.write(
                        " {:2s} {:6d}{:14.6e}{:14.6e}{:14.6e}{:14.6e}{:14.6e}\n".format(
                            data.elem[iel].bcs[ibc, iface][0],
                            data.elem[iel].bcs[ibc, iface][1],
                            data.elem[iel].bcs[ibc, iface][3],
                            data.elem[iel].bcs[ibc, iface][4],
                            data.elem[iel].bcs[ibc, iface][5],
                            data.elem[iel].bcs[ibc, iface][6],
                            data.elem[iel].bcs[ibc, iface][7],
                        )
                    )
                else:
                    outfile.write(
                        " {:2s} {:11d}{:1d}{:18.11e}{:18.11e}{:18.11e}{:18.11e}{:18.11e}\n".format(
                            data.elem[iel].bcs[ibc, iface][0],
                            data.elem[iel].bcs[ibc, iface][1],
                            data.elem[iel].bcs[ibc, iface][2],
                            data.elem[iel].bcs[ibc, iface][3],
                            data.elem[iel].bcs[ibc, iface][4],
                            data.elem[iel].bcs[ibc, iface][5],
                            data.elem[iel].bcs[ibc, iface][6],
                            data.elem[iel].bcs[ibc, iface][7],
                        )
                    )

    if data.nbc < 2:
        outfile.write("  ***** NO THERMAL BOUNDARY CONDITIONS *****\n")
    outfile.write("    0 PRESOLVE/RESTART OPTIONS  *****\n")
    outfile.write("    7         INITIAL CONDITIONS *****\n")
    outfile.write(" C Default\n")
    outfile.write(" C Default\n")
    outfile.write(" C Default\n")
    outfile.write(" C Default\n")
    outfile.write(" C Default\n")
    outfile.write(" C Default\n")
    outfile.write(" C Default\n")
    outfile.write("  ***** DRIVE FORCE DATA ***** BODY FORCE, FLOW, Q\n")
    outfile.write("    4                 Lines of Drive force data follow\n")
    outfile.write(" C\n")
    outfile.write(" C\n")
    outfile.write(" C\n")
    outfile.write(" C\n")
    outfile.write(" ***** Variable Property Data ***** Overrrides Parameter data.\n")
    outfile.write("    1 Lines follow.\n")
    outfile.write("    0 PACKETS OF DATA FOLLOW\n")
    outfile.write(" ***** HISTORY AND INTEGRAL DATA *****\n")
    outfile.write("    0   POINTS.  Hcode, I,J,H,IEL\n")
    outfile.write(" ***** OUTPUT FIELD SPECIFICATION *****\n")
    outfile.write("    6 SPECIFICATIONS FOLLOW\n")
    outfile.write("    F      COORDINATES\n")
    outfile.write("    T      VELOCITY\n")
    outfile.write("    T      PRESSURE\n")
    outfile.write("    F      TEMPERATURE\n")
    outfile.write("    F      TEMPERATURE GRADIENT\n")
    outfile.write("    0      PASSIVE SCALARS\n")
    outfile.write(" ***** OBJECT SPECIFICATION *****\n")
    outfile.write("        0 Surface Objects\n")
    outfile.write("        0 Volume  Objects\n")
    outfile.write("        0 Edge    Objects\n")
    outfile.write("        0 Point   Objects\n")
    #
    # close file
    outfile.close()
    #
    # output
    return 0


def readre2(fname):
    """A function for reading .re2 files for nek5000

    Parameters
    ----------
    fname : str
            file name
    """
    #
    infile = open(fname, "rb")
    # the header for re2 files is 80 ASCII bytes, something like
    # #v002    18669  2    18669 this is the hdr                                      %
    header = infile.read(80).split()
    nel = int(header[1])
    ndim = int(header[2])
    # always double precision
    wdsz = 8
    realtype = "d"

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
        raise ValueError("Could not interpret endianness")

    #
    # there are no GLL points here, only quad/hex vertices
    lr1 = [2, 2, ndim - 1]
    npel = 2**ndim
    # the file only contains geometry
    var = [ndim, 0, 0, 0, 0]
    # allocate structure
    data = HexaData(ndim, nel, lr1, var, 1)
    #
    # some metadata
    data.wdsz = wdsz
    data.endian = endian
    #
    # read the whole geometry into a buffer
    # for some reason each element is prefixed with 8 bytes of zeros, it's not clear to me why.
    # This is the reason for the +1 here, then the first number is ignored.
    buf = infile.read((ndim * npel + 1) * wdsz * nel)
    # elem_shape = [ndim, ndim-1, 2, 2]  # nvar, lz, ly, lx
    for iel, el in enumerate(data.elem):
        fi = np.frombuffer(
            buf,
            dtype=emode + realtype,
            count=ndim * npel + 1,
            offset=(ndim * npel + 1) * wdsz * iel,
        )
        # the data is stored in the following order (2D|3D):
        # x1, x2, x4, x3; | x5, x6, x8, x7;
        # y1, y2, y4, y3; | y5, y6, y8, y7;
        # ----------------
        # z1, z2, z4, z3; z5, z6, z8, z7;
        # where 1-8 is the ordering of the points in memory
        for idim in range(ndim):  # x, y, [z]
            for iz in range(
                ndim - 1
            ):  # this does only one iteration in 2D, and in 3D does one iteration for 1-4 and one for 5-8
                el.pos[idim, iz, 0, 0] = fi[npel * idim + 4 * iz + 1]
                el.pos[idim, iz, 0, 1] = fi[npel * idim + 4 * iz + 2]
                el.pos[idim, iz, 1, 1] = fi[npel * idim + 4 * iz + 3]
                el.pos[idim, iz, 1, 0] = fi[npel * idim + 4 * iz + 4]
    #
    # read curved sides
    # the number of curved sides is stored as a double,
    # then each curved side is 64 bytes:
    # iel iface p1 p2 p3 p4 p5 ctype
    # where p1-5 are f64 parameters and ctype is the type of curvature in ASCII
    ncparam = 8
    buf = infile.read(wdsz)
    ncurv = int(np.frombuffer(buf)[0])
    logger.debug(f"Found {ncurv} curved sides")
    data.ncurv = ncurv
    buf = infile.read(wdsz * (ncparam * ncurv))
    for icurv in range(ncurv):
        # interpret the data
        curv = np.frombuffer(
            buf,
            dtype=emode + realtype,
            count=ncparam,
            offset=icurv * ncparam * wdsz,
        )
        iel = int(curv[0]) - 1
        iedge = int(curv[1]) - 1
        cparams = curv[2:7]
        # select only the first byte, because it turns out the later bytes may contain garbage.
        # typically, it's b'C\x00\x00\x00\x00\x00\x00\x00' or b'C\x00\x00\x00\x00\x00\xe0?'.
        # AFAIK the curvature types are always one character long anyway.
        ctype = curv[7].tobytes()[:1].decode("utf-8")
        # fill in the data structure
        data.elem[iel].curv[iedge, :] = cparams
        data.elem[iel].ccurv[iedge] = ctype
    #
    # read boundary conditions
    # there can be more than one field, and we won't know until we'vre reached the end
    nbcparam = 8
    buf = infile.read(wdsz)
    ifield = 0
    while buf != b"":
        # the data is initialized with one BC field, we might need to allocate another
        if ifield > 0:
            data.nbc = data.nbc + 1
            for el in data.elem:
                empty_bcs = np.zeros(el.bcs[:1, :].shape, dtype=el.bcs.dtype)
                el.bcs = np.concatenate((el.bcs, empty_bcs))
        nbclines = int(np.frombuffer(buf)[0])
        logger.debug(
            f"Found {nbclines} external boundary conditions for field {ifield}"
        )
        buf = infile.read(wdsz * (nbcparam * nbclines))
        for ibc in range(nbclines):
            # interpret the data
            bc = np.frombuffer(
                buf,
                dtype=emode + realtype,
                count=nbcparam,
                offset=ibc * nbcparam * wdsz,
            )
            iel = int(bc[0]) - 1
            iface = int(bc[1]) - 1
            bcparams = bc[2:7]
            bctype = bc[7].tobytes().decode("utf-8").rstrip()  # remove trailing spaces
            # fill in the data structure
            data.elem[iel].bcs[ifield, iface][0] = bctype
            data.elem[iel].bcs[ifield, iface][1] = iel + 1
            data.elem[iel].bcs[ifield, iface][2] = iface + 1
            for ipar in range(5):
                data.elem[iel].bcs[ifield, iface][3 + ipar] = bcparams[ipar]
        ifield = ifield + 1
        # try reading the number of conditions in the next field
        buf = infile.read(wdsz)
    infile.close()
    return data


def writere2(fname, data):
    """A function for writing binary .re2 files for nek5000

    Parameters
    ----------
    fname : str
            file name
    data : :class:`pymech.core.HexaData`
            data structure
    """
    #
    # ---------------------------------------------------------------------------
    # CHECK INPUT DATA
    # ---------------------------------------------------------------------------
    #
    # We could extract the corners, but for now just return an error if lr1 is too large
    if data.lr1 != [2, 2, data.ndim - 1]:
        raise ValueError(
            "wrong element dimensions for re2 file! {} != {}".format(
                data.lr1, [2, 2, data.ndim - 1]
            )
        )
    #
    if data.var[0] != data.ndim:
        raise ValueError(
            "wrong number of geometric variables for re2 file! expected {}, found {}".format(
                data.ndim, data.var[0]
            )
        )
    #
    # Open file
    outfile = open(fname, "wb")
    #
    # ---------------------------------------------------------------------------
    # WRITE HEADER
    # ---------------------------------------------------------------------------
    #
    # always double precision
    #  wdsz = 8
    #  realtype = 'd'
    nel = data.nel
    ndim = data.ndim
    header = f"#v002{nel:9d}{ndim:3d}{nel:9d} this is the hdr"
    header = header.ljust(80)
    outfile.write(header.encode("utf-8"))
    #
    # decide endianness
    if data.endian in ("big", "little"):
        byteswap = data.endian != sys.byteorder
        logger.debug(f"Writing {data.endian}-endian file")
    else:
        byteswap = False
        logger.warning(
            f"Unrecognized endianness {data.endian}, "
            f"writing native {sys.byteorder}-endian file"
        )

    def correct_endianness(a):
        """Return the array with the requested endianness"""
        if byteswap:
            return a.byteswap()
        else:
            return a

    #
    # write tag (to specify endianness)
    endianbytes = np.array([6.54321], dtype=np.float32)
    correct_endianness(endianbytes).tofile(outfile)
    #
    # ---------------------------------------------------------------------------
    # WRITE DATA
    # ---------------------------------------------------------------------------
    #
    # compute total number of points per element
    npel = 2**ndim

    def write_data_to_file(a):
        """Write the geometry of an element to the output file in double precision"""
        correct_endianness(a).tofile(outfile)

    #
    # write geometry (adding eight bytes of zeros before each element)
    xyz = np.zeros(
        (npel * ndim + 1,)
    )  # array storing reordered geometry data (with a zero in the first position)
    for el in data.elem:
        # the data is stored in the following order (2D|3D):
        # x1, x2, x4, x3; | x5, x6, x8, x7;
        # y1, y2, y4, y3; | y5, y6, y8, y7;
        # ----------------
        # z1, z2, z4, z3; z5, z6, z8, z7;
        # where 1-8 is the ordering of the points in memory
        for idim in range(ndim):  # x, y, [z]
            for iz in range(
                ndim - 1
            ):  # this does only one iteration in 2D, and in 3D does one iteration for 1-4 and one for 5-8
                xyz[npel * idim + 4 * iz + 1] = el.pos[idim, iz, 0, 0]
                xyz[npel * idim + 4 * iz + 2] = el.pos[idim, iz, 0, 1]
                xyz[npel * idim + 4 * iz + 3] = el.pos[idim, iz, 1, 1]
                xyz[npel * idim + 4 * iz + 4] = el.pos[idim, iz, 1, 0]
        write_data_to_file(xyz)
    #
    # write curve sides data
    # locate curved edges
    curved_edges = []
    for iel, el in enumerate(data.elem):
        for iedge in range(12):
            if el.ccurv[iedge] != "":
                curved_edges.append((iel, iedge))
    # write number of curved edges
    ncurv = len(curved_edges)
    if ncurv != data.ncurv:
        logger.warning(
            f"wrong number of curved edges: expected {data.ncurv}, found {ncurv}"
        )
    ncurvf = np.array([ncurv], dtype=np.float64)
    write_data_to_file(ncurvf)
    # format curve data
    cdata = np.zeros((ncurv,), dtype="f8, f8, f8, f8, f8, f8, f8, S8")
    for cdat, (iel, iedge) in zip(cdata, curved_edges):
        el = data.elem[iel]
        cdat[0] = iel + 1
        cdat[1] = iedge + 1
        # curve parameters
        for j in range(5):
            cdat[2 + j] = el.curv[iedge, j]
        # encode the string as a byte array padded with spaces
        cdat[7] = el.ccurv[iedge].encode("utf-8")
    # write to file
    write_data_to_file(cdata)
    #
    # write boundary conditions for each field
    for ifield in range(data.nbc):
        # locate faces with boundary conditions
        bc_faces = []
        for iel, el in enumerate(data.elem):
            for iface in range(2 * ndim):
                bctype = el.bcs[ifield, iface][0]
                # internal boundary conditions are not written to .re2 files by reatore2
                # and are apparently ignored by Nek5000 even in .rea files
                if bctype != "" and bctype != "E":
                    bc_faces.append((iel, iface))
        nbcs = len(bc_faces)
        nbcsf = np.array([nbcs], dtype=np.float64)
        write_data_to_file(nbcsf)
        # initialize and format data
        bcdata = np.zeros((nbcs,), dtype="f8, f8, f8, f8, f8, f8, f8, S8")
        for bc, (iel, iface) in zip(bcdata, bc_faces):
            el = data.elem[iel]
            bc[0] = iel + 1
            bc[1] = iface + 1
            for j in range(5):
                bc[2 + j] = el.bcs[ifield, iface][3 + j]
            # encode the string as a byte array padded with spaces
            bc[7] = el.bcs[ifield, iface][0].encode("utf-8").ljust(8)
        # write to file
        write_data_to_file(bcdata)
    #
    # close file
    outfile.close()
    # rerurn
    return 0
