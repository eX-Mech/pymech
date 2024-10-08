"""Core data structures for pymech"""

import copy
import itertools
from functools import partial, reduce
from itertools import product
from textwrap import dedent, indent

import numpy as np

from pymech.log import logger

"""Repeat N times. Pythonic idiom to use when the iterated value is discarded.

Example
-------
Instead of:

>>> [0 for _ in range(10)]

You could use:

>>> [0 for _ in repeat(10)]

"""
repeat = partial(itertools.repeat, None)


# ==============================================================================
class DataLims:
    """A class containing the extrema of all quantities stored in the mesh

    Attributes
    ----------
    - pos:  x,y,z   min,max
    - vel:  u,v,w   min,max
    - pres: p       min,max
    - temp: T       min,max
    - scal: s_i     min,max

    """

    def __init__(self, elements):
        self._variables = ("pos", "vel", "pres", "temp", "scal")

        aggregated_lims = reduce(self._lims_aggregator, elements)
        for var in self._variables:
            agg_lims_var = aggregated_lims[var]
            # set minimum, maximum of variables as a nested tuple
            setattr(
                self,
                var,
                tuple((float(lims[0]), float(lims[1])) for lims in zip(*agg_lims_var)),
            )

        # prevent further mutation of attributes via __setattr__
        self._initialized = True

    def __repr__(self):
        return dedent(
            f"""\
          * x:         {self.pos[0]}
          * y:         {self.pos[1]}
          * z:         {self.pos[2]}"""
        )

    def __setattr__(self, name, value):
        if hasattr(self, "_initialized") and self._initialized:
            raise AttributeError(f"Setting attribute {name} is not permitted")
        else:
            super().__setattr__(name, value)

    def _lims_per_element(self, elem):
        """Get local limits for a given element."""
        if isinstance(elem, dict):
            return elem

        axis = (1, 2, 3)
        elem_lims = {
            var: (getattr(elem, var).min(axis), getattr(elem, var).max(axis))
            for var in self._variables
        }
        return elem_lims

    def _lims_aggregator(self, elem1, elem2):
        """Reduce local limits to global limits."""
        l1 = self._lims_per_element(elem1)
        l2 = self._lims_per_element(elem2)

        aggregated_lims = {
            var: (
                np.minimum(l1[var][0], l2[var][0]),
                np.maximum(l1[var][1], l2[var][1]),
            )
            for var in self._variables
        }
        return aggregated_lims


# ==============================================================================
class Elem:
    """A class containing one hexahedral element of Nek5000/SIMSON flow
    field.

    Parameters
    ----------
    var : iterable
        Iterable of integers of size 5, indicating how many variables are to be initialized
    lr1 : iterable
        Iterable of integers of size 3, defining the shape of an element as ``(lx, ly, lz)``
    nbc : int
        Number of boundary conditions
    dtype : str
        Floating point data type. Typical values are 'f4' or 'float32' for
        single precision, 'f8' or 'float64' for double precision

    """

    def __init__(self, var, lr1, nbc, dtype="float64"):
        #                    x,y,z   lz      ly      lx
        self.pos = np.zeros((3, lr1[2], lr1[1], lr1[0]), dtype=dtype)
        #                    one per edge
        self.curv = np.zeros((12, 5), dtype=dtype)
        #             curvature type
        self.ccurv = ["" for _ in repeat(12)]
        #                    u,v,w   lz      ly      lx
        self.vel = np.zeros((3, lr1[2], lr1[1], lr1[0]), dtype=dtype)
        #                    p       lz      ly      lx
        self.pres = np.zeros((var[2], lr1[2], lr1[1], lr1[0]), dtype=dtype)
        #                    T       lz      ly      lx
        self.temp = np.zeros((var[3], lr1[2], lr1[1], lr1[0]), dtype=dtype)
        #                    s_i     lz      ly      lx
        self.scal = np.zeros((var[4], lr1[2], lr1[1], lr1[0]), dtype=dtype)
        #                    list of 8 parameters, one per face
        #                    one column for velocity, one for temperature, and one for each scalar
        self.bcs = np.zeros((nbc, 6), dtype="U3, i4, i4" + f", {dtype}" * 5)

    def __repr__(self):
        message = f"<elem centered at {self.centroid}>"
        return message

    @property
    def centroid(self):
        return self.pos.mean(axis=(1, 2, 3))

    def smallest_edge(self):
        """returns the length of the smallest edge, neglecting curvature"""

        # get coordinates of points
        x1, y1, z1 = self.pos[:, 0, 0, 0]
        x2, y2, z2 = self.pos[:, 0, 0, 1]
        x3, y3, z3 = self.pos[:, 0, 1, 0]
        x4, y4, z4 = self.pos[:, 0, 1, 1]
        # compute squares of edges lengths
        edges_l2 = np.zeros((12,))
        edges_l2[0] = (x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2
        edges_l2[1] = (x3 - x2) ** 2 + (y3 - y2) ** 2 + (z3 - z2) ** 2
        edges_l2[2] = (x4 - x3) ** 2 + (y4 - y3) ** 2 + (z4 - z3) ** 2
        edges_l2[3] = (x1 - x4) ** 2 + (y1 - y4) ** 2 + (z1 - z4) ** 2
        # the dimension is not stored but we can cheat
        if self.pos.shape[1] > 1:
            ndim = 3
        else:
            ndim = 2
        if ndim > 2:
            # in 3D, do the same for the upper face, and also the side edges
            x5, y5, z5 = self.pos[:, 1, 0, 0]
            x6, y6, z6 = self.pos[:, 1, 0, 1]
            x7, y7, z7 = self.pos[:, 1, 1, 0]
            x8, y8, z8 = self.pos[:, 1, 1, 1]
            edges_l2[4] = (x6 - x5) ** 2 + (y6 - y5) ** 2 + (z6 - z5) ** 2
            edges_l2[5] = (x7 - x6) ** 2 + (y7 - y6) ** 2 + (z7 - z6) ** 2
            edges_l2[6] = (x8 - x7) ** 2 + (y8 - y7) ** 2 + (z8 - z7) ** 2
            edges_l2[7] = (x5 - x8) ** 2 + (y5 - y8) ** 2 + (z5 - z8) ** 2
            edges_l2[8] = (x5 - x1) ** 2 + (y5 - y1) ** 2 + (z5 - z1) ** 2
            edges_l2[9] = (x6 - x2) ** 2 + (y6 - y2) ** 2 + (z6 - z2) ** 2
            edges_l2[10] = (x7 - x3) ** 2 + (y7 - y3) ** 2 + (z7 - z3) ** 2
            edges_l2[11] = (x8 - x4) ** 2 + (y8 - y4) ** 2 + (z8 - z4) ** 2
            return np.sqrt(edges_l2.min())
        else:
            return np.sqrt(edges_l2[:4].min())

    def face_center(self, i):
        """Return the coordinates (x, y, z) of the center of the face number i"""

        if i == 0:
            kx1, ky1, kz1 = 0, 0, 0
            kx2, ky2, kz2 = 0, 0, -1
            kx3, ky3, kz3 = -1, 0, 0
            kx4, ky4, kz4 = -1, 0, -1
        elif i == 1:
            kx1, ky1, kz1 = -1, 0, 0
            kx2, ky2, kz2 = -1, 0, -1
            kx3, ky3, kz3 = -1, -1, 0
            kx4, ky4, kz4 = -1, -1, -1
        elif i == 2:
            kx1, ky1, kz1 = 0, -1, 0
            kx2, ky2, kz2 = 0, -1, -1
            kx3, ky3, kz3 = -1, -1, 0
            kx4, ky4, kz4 = -1, -1, -1
        elif i == 3:
            kx1, ky1, kz1 = 0, 0, 0
            kx2, ky2, kz2 = 0, 0, -1
            kx3, ky3, kz3 = 0, -1, 0
            kx4, ky4, kz4 = 0, -1, -1
        elif i == 4:
            kx1, ky1, kz1 = 0, 0, 0
            kx2, ky2, kz2 = 0, -1, 0
            kx3, ky3, kz3 = -1, 0, 0
            kx4, ky4, kz4 = -1, -1, 0
        elif i == 5:
            kx1, ky1, kz1 = 0, 0, -1
            kx2, ky2, kz2 = 0, -1, -1
            kx3, ky3, kz3 = -1, 0, -1
            kx4, ky4, kz4 = -1, -1, -1
        else:
            logger.error(f"Invalid face number {i} (must be between 0 and 5)")
        (x1, y1, z1) = self.pos[:, kz1, ky1, kx1]
        (x2, y2, z2) = self.pos[:, kz2, ky2, kx2]
        (x3, y3, z3) = self.pos[:, kz3, ky3, kx3]
        (x4, y4, z4) = self.pos[:, kz4, ky4, kx4]
        return (
            0.25 * (x1 + x2 + x3 + x4),
            0.25 * (y1 + y2 + y3 + y4),
            0.25 * (z1 + z2 + z3 + z4),
        )


# ==============================================================================
class HexaData:
    """A class containing data related to a hexahedral mesh"""

    def __init__(self, ndim, nel, lr1, var, nbc=0, dtype="float64"):
        self.ndim = ndim
        self.nel = nel
        self.ncurv = []
        self.nbc = nbc
        self.var = var
        self.lr1 = lr1
        self.time = []
        self.istep = []
        self.wdsz = []
        self.endian = []
        if isinstance(dtype, type):
            # For example np.float64 -> "float64"
            dtype = dtype.__name__

        self.elem = [Elem(var, lr1, nbc, dtype) for _ in repeat(nel)]
        self.elmap = np.linspace(1, nel, nel, dtype=np.int32)

    def __repr__(self):
        representation = dedent(
            f"""\
        <pymech.core.HexaData>
        Dimensions:    {self.ndim}
        Precision:     {self.wdsz} bytes
        Mesh limits:\n{indent(repr(self.lims), " "*10)}
        Time:
          * time:      {self.time}
          * istep:     {self.istep}
        Elements:
          * nel:       {self.nel}
          * elem:      [{self.elem[0]}
                        ...
                        {self.elem[-1]}]
        """
        )

        return representation

    @property
    def lims(self):
        return DataLims(self.elem)

    def check_connectivity(self, tol=1e-3):
        """Check element connectivity, specifically for matching boundary
        conditions and geometry. Errors are reported as logging messages.

        Parameters
        ----------
        tol : float
            relative tolerance (compared to the smallest edge of adjacent elements)
            for detecting whether faces are at the same location

        """
        dim = self.ndim
        err = False
        for (iel, el), ibc, iface in product(
            enumerate(self.elem), range(self.nbc), range(2 * dim)
        ):
            cbc = el.bcs[ibc, iface][0]
            if cbc == "E" or cbc == "P":
                connected_iel = int(el.bcs[ibc, iface][3]) - 1
                connected_face = int(el.bcs[ibc, iface][4]) - 1
                xc, yc, zc = el.face_center(iface)
                if connected_iel < 0 or connected_iel >= self.nel:
                    err = True
                    logger.error(
                        f"face {iface} of element {iel} is connected ('{cbc}') to face "
                        f"{connected_face} of the nonexistent element {connected_iel}"
                    )
                    logger.error(f"face center: ({xc:.6e} {yc:.6e} {zc:.6e})")
                else:
                    cbc1 = self.elem[connected_iel].bcs[ibc, connected_face][0]
                    iel1 = int(self.elem[connected_iel].bcs[ibc, connected_face][3]) - 1
                    iface1 = (
                        int(self.elem[connected_iel].bcs[ibc, connected_face][4]) - 1
                    )
                    xc1, yc1, zc1 = self.elem[connected_iel].face_center(connected_face)
                    if cbc1 != cbc or iel1 != iel or iface1 != iface:
                        err = True
                        logger.error(
                            "mismatched boundary conditions: "
                            f"face {iface + 1} of element {iel + 1} with "
                            f"condition '{cbc}' is connected to face {connected_face + 1} "
                            f"of element {connected_iel + 1}, which has condition '{cbc1}' "
                            f"and points to face {iface1} of element {iel1}"
                        )
                        logger.error(
                            f"face centers: ({xc:.6e} {yc:.6e} {zc:.6e}), ({xc1:.6e} {yc1:.6e} {zc1:.6e})"
                        )
                    elif cbc == "E":  # no check for 'P' yet, but it should be possible
                        dist = np.sqrt(
                            (xc - xc1) ** 2 + (yc - yc1) ** 2 + (zc - zc1) ** 2
                        )
                        max_dist = tol * min(
                            el.smallest_edge(), self.elem[connected_iel].smallest_edge()
                        )
                        if dist > max_dist:
                            err = True
                            logger.error(
                                "mismatched face locations: "
                                f"face {iface + 1} of element {iel + 1} and "
                                f"face {connected_face + 1} of element {connected_iel + 1} "
                                "are connected but are not at the same location"
                            )
                            logger.error(
                                f"face centers: ({xc:.6e} {yc:.6e} {zc:.6e}), ({xc1:.6e} {yc1:.6e} {zc1:.6e})"
                            )

        if err:
            raise ValueError(
                "Some errors were encountered while checking connectivity."
            )

        return not err

    def check_bcs_present(self):
        """
        Returns True if and only if all faces of all elements have boundary conditions applied.

        Note that this function returning False does not mean the mesh is invalid: it is not mandatory
        to define internal boundary conditions for Nek5000.
        """

        res = True
        for (iel, el), ibc, iface in product(
            enumerate(self.elem), range(self.nbc), range(2 * self.ndim)
        ):
            if el.bcs[ibc, iface][0] == "":
                res = False
                logger.error(
                    f"missing boundary condition at element {iel}, face {iface}, field {ibc}"
                )

        if not res:
            raise ValueError(
                "Some errors were encountered while checking boundary conditions."
            )

        return res

    def merge(self, other, tol=1e-2, ignore_empty=True, ignore_all_bcs=False):
        """
        Merges another :class:`pymech.core.HexaData` into the current one and connects it

        Parameters
        ----------
        other: :class:`pymech.core.HexaData`
                mesh to merge into self

        tol: float
                maximum distance, relative to the smallest edge of neighbouring elements, at which faces are considered touching

        ignore_empty: bool
                if True, the faces with an empty boundary condition ('') will be treated as internal faces and will not be merged.
                This is useful if internal boundary conditions are not defined and will make the operation *much* faster,
                but requires boundary conditions to be defined on the faces to be merged.

        ignore_all_bcs: bool
                if True, the boundary conditions will not be changed.
                This is likely to result in invalid boundary conditions at the interface between the merged meshes.
                This option is intended for fast merging in a situation in which the boundary conditions
                are either irrelevant or will be defined or corrected later.
        """

        # perform some consistency checks
        if self.ndim != other.ndim:
            raise ValueError(
                f"Cannot merge meshes of dimensions {self.ndim} and {other.ndim}!"
            )
        if self.lr1[0] != other.lr1[0]:
            raise ValueError(
                "Cannot merge meshes of different polynomial orders ({} != {})".format(
                    self.lr1[0], other.lr1[0]
                )
            )

        # add the new elements (in an inconsistent state if there are internal boundary conditions)
        nel1 = self.nel
        self.nel = self.nel + other.nel
        self.ncurv = self.ncurv + other.ncurv
        self.elmap = np.concatenate((self.elmap, other.elmap + nel1))
        # the deep copy is required here to avoid leaving the 'other' mesh in an inconsistent state by modifying its elements
        mesh_add = copy.deepcopy(other)
        mesh_add.offset_connectivity(nel1)
        self.elem = self.elem + mesh_add.elem

        # check how many boundary condition fields we have
        nbc = min(self.nbc, other.nbc)

        # glue common faces together
        # only look for the neighbour in the first BC field because it should be the same in all fields.
        nfaces = 2 * self.ndim
        nchanges = 0  # counter for the boundary conditions connected
        if nbc == 0 or ignore_all_bcs:
            # Quickly exit the function
            logger.warning("No pairs of faces to merge.")
            return nchanges

        for iel0, iface0 in product(range(nel1, self.nel), range(nfaces)):
            elem0 = self.elem[iel0]
            bc0 = elem0.bcs[0, iface0][0]

            if bc0 != "E" and not (ignore_empty and bc0 == ""):
                # boundary element, look if it can be connected to something
                for iel1, iface1 in product(range(nel1), range(nfaces)):
                    elem1 = self.elem[iel1]
                    bc1 = elem1.bcs[0, iface1][0]

                    if bc1 != "E" and not (ignore_empty and bc1 == ""):
                        # if the centers of the faces are close, connect them together
                        x0, y0, z0 = elem0.face_center(iface0)
                        x1, y1, z1 = elem1.face_center(iface1)
                        dist = np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2 + (z1 - z0) ** 2)
                        dist_ref = min(elem0.smallest_edge(), elem1.smallest_edge())
                        if dist <= tol * dist_ref:
                            # reconnect the periodic faces together (assumes that all fields are periodic)
                            if bc0 == bc1 == "P":
                                iel_p0 = int(elem0.bcs[0, iface0][3]) - 1
                                iel_p1 = int(elem1.bcs[0, iface1][3]) - 1
                                iface_p0 = int(elem0.bcs[0, iface0][4]) - 1
                                iface_p1 = int(elem1.bcs[0, iface1][4]) - 1
                                for ibc in range(nbc):
                                    elem_p0_bcs = self.elem[iel_p0].bcs[ibc, iface_p0]
                                    elem_p1_bcs = self.elem[iel_p1].bcs[ibc, iface_p1]

                                    elem_p0_bcs[0] = "P"
                                    elem_p1_bcs[0] = "P"
                                    elem_p0_bcs[3] = iel_p1 + 1
                                    elem_p1_bcs[3] = iel_p0 + 1
                                    elem_p0_bcs[4] = iface_p1 + 1
                                    elem_p1_bcs[4] = iface_p0 + 1

                            for ibc in range(nbc):
                                elem0.bcs[ibc, iface0][0] = "E"
                                elem1.bcs[ibc, iface1][0] = "E"
                                elem0.bcs[ibc, iface0][3] = iel1 + 1
                                elem1.bcs[ibc, iface1][3] = iel0 + 1
                                elem0.bcs[ibc, iface0][4] = iface1 + 1
                                elem1.bcs[ibc, iface1][4] = iface0 + 1
                            nchanges = nchanges + 1

        logger.debug(f"merged {nchanges} pairs of faces")
        return nchanges

    def get_points(self):
        """
        Returns an array containing the coordinates of all the points in the mesh as a (nel, lx1*ly1*lz1, 3) array
        """

        lx1, ly1, lz1 = self.lr1
        nxyz = lx1 * ly1 * lz1
        xyz = np.zeros((self.nel, nxyz, 3))
        for el, lxyz in zip(self.elem, xyz):
            for idim in range(3):
                lxyz[:, idim] = el.pos[idim, ...].ravel()

        return xyz

    def update_ncurv(self):
        """
        Updates the metadata `ncurv` integer to match the actual number of curved faces present in the mesh
        """

        ncurv = 0
        for el in self.elem:
            for iedge in range(12):
                if el.ccurv[iedge] != "":
                    ncurv = ncurv + 1
        self.ncurv = ncurv

    def offset_connectivity(self, offset: int, iel_min=0):
        """
        Adds a value to the index of the elements connected via internal or periodic
        boundary conditions to elements of the mesh. This is used to keep the connectivity
        valid when deleting or inserting elements in the mesh.

        Parameters
        ----------
        offset  : int
            The value by which to offset the indices
        iel_min : int
            The first element (in zero-based indexing) to offset
        """

        for el, ibc, iface in product(self.elem, range(self.nbc), range(2 * self.ndim)):
            bc = el.bcs[ibc, iface][0]
            if bc == "E" or bc == "P":
                if (
                    int(el.bcs[ibc, iface][3]) > iel_min
                ):  # the connected element number is 1-indexed
                    el.bcs[ibc, iface][3] += offset

        # also fix the index of the elements themselves in its their BC if relevant
        for el, ibc, iface in product(
            self.elem[iel_min:], range(self.nbc), range(2 * self.ndim)
        ):
            if int(el.bcs[ibc, iface][1]) > iel_min:
                el.bcs[ibc, iface][1] += offset
