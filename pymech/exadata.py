"""Data structures for pymech"""
import copy
from textwrap import dedent, indent
from itertools import product
from functools import reduce

import numpy as np
from pymech.log import logger


# ==============================================================================
class datalims:
    """A class containing the extrema of all quantities stored in the mesh"""

    def __init__(self, nb_var, elements):
        #                    x,y,z   min,max
        self.pos = np.zeros((3, 2))
        #                    u,v,w   min,max
        self.vel = np.zeros((3, 2))
        #                    p       min,max
        self.pres = np.zeros((nb_var[2], 2))
        #                    T       min,max
        self.temp = np.zeros((nb_var[3], 2))
        #                    s_i     min,max
        self.scal = np.zeros((nb_var[4], 2))
        #
        self._variables = ("pos", "vel", "pres", "temp", "scal")

        aggregated_lims = reduce(self._lims_aggregator, elements)
        for var in self._variables:
            agg_lims_var = aggregated_lims[var]
            # set minimum
            getattr(self, var)[:, 0] = agg_lims_var[0]
            # set maximum
            getattr(self, var)[:, 1] = agg_lims_var[1]

    def __repr__(self):
        return dedent(
            f"""\
          * x:         {self.pos[0]}
          * y:         {self.pos[1]}
          * z:         {self.pos[2]}"""
        )

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
class elem:
    """A class containing one nek element/SIMSON flow field"""

    def __init__(self, var, lr1, nbc):
        #                    x,y,z   lz      ly      lx
        self.pos = np.zeros((3, lr1[2], lr1[1], lr1[0]))
        #                    one per edge
        self.curv = np.zeros((12, 5))
        #             curvature type
        self.ccurv = ["" for i in range(12)]
        #                    u,v,w   lz      ly      lx
        self.vel = np.zeros((3, lr1[2], lr1[1], lr1[0]))
        #                    p       lz      ly      lx
        self.pres = np.zeros((var[2], lr1[2], lr1[1], lr1[0]))
        #                    T       lz      ly      lx
        self.temp = np.zeros((var[3], lr1[2], lr1[1], lr1[0]))
        #                    s_i     lz      ly      lx
        self.scal = np.zeros((var[4], lr1[2], lr1[1], lr1[0]))
        #                    list of 8 parameters, one per face
        #                    one column for velocity, one for temperature, and one for each scalar
        self.bcs = np.zeros((nbc, 6), dtype="U3, i4, i4, f8, f8, f8, f8, f8")

    def __repr__(self):
        message = f"<elem centered at {self.centroid}>"
        return message

    @property
    def centroid(self):
        return self.pos.mean(axis=(1, 2, 3))

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
class exadata:
    """A class containing data for reading/writing binary simulation files"""

    def __init__(self, ndim, nel, lr1, var, nbc=0):
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
        self.elem = [elem(var, lr1, nbc) for i in range(nel)]
        self.elmap = np.linspace(1, nel, nel, dtype=np.int32)

    def __repr__(self):
        representation = dedent(
            f"""\
        <pymech.exadata.exadata>
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
        return datalims(self.var, self.elem)

    def check_connectivity(self):
        """Check element connectivity, specifically for matching boundary
        conditions and geometry. Errors are reported as logging messages.

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
                cbc1 = self.elem[connected_iel].bcs[ibc, connected_face][0]
                iel1 = int(self.elem[connected_iel].bcs[ibc, connected_face][3]) - 1
                iface1 = int(self.elem[connected_iel].bcs[ibc, connected_face][4]) - 1
                if iel1 < 0 or iel1 >= self.nel:
                    err = True
                    logger.error(
                        f"face {iface} of element {iel} is connected to face "
                        f"{connected_face} of the nonexistent element {connected_iel}"
                    )
                else:
                    if cbc1 != cbc:
                        err = True
                        logger.error(
                            "mismatched boundary conditions: "
                            f"face {iface + 1} of element {iel + 1} with "
                            f"condition {cbc} is connected to face {connected_face + 1} "
                            f"of element {connected_iel + 1} with condition {cbc1}"
                        )
                    if iel1 != iel:
                        err = True
                        logger.error(
                            "mismatched elements: "
                            f"face {iface + 1} of element {iel + 1} "
                            f"is connected to face {connected_face + 1} "
                            f"of element {connected_iel + 1} but that face is "
                            f"connected to face {iface + 1} of element {iel + 1}"
                        )
                    if iface1 != iface:
                        err = True
                        logger.error(
                            "mismatched faces: "
                            f"face {iface + 1} of element {iel + 1} "
                            f"is connected to face {connected_face + 1} "
                            f"of element {connected_iel + 1} but that face is "
                            f"connected to face {iface + 1} of element {iel + 1}"
                        )
        return not err

    def merge(self, other, tol=1e-9, ignore_empty=True, ignore_all_bcs=False):
        """
        Merges another exadata into the current one and connects it

        Parameters
        ----------
        other: exadata
                mesh to merge into self

        tol: float
                maximum distance at which points are considered identical

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
            logger.error(
                f"Cannot merge meshes of dimensions {self.ndim} and {other.ndim}!"
            )
            return -1
        if self.lr1[0] != other.lr1[0]:
            logger.error(
                "Cannot merge meshes of different polynomial orders ({} != {})".format(
                    self.lr1[0], other.lr1[0]
                )
            )
            return -2

        # add the new elements (in an inconsistent state if there are internal boundary conditions)
        nel1 = self.nel
        self.nel = self.nel + other.nel
        self.ncurv = self.ncurv + other.ncurv
        self.elmap = np.concatenate((self.elmap, other.elmap + nel1))
        # the deep copy is required here to avoid leaving the 'other' mesh in an inconsistent state by modifying its elements
        self.elem = self.elem + copy.deepcopy(other.elem)

        # check how many boundary condition fields we have
        nbc = min(self.nbc, other.nbc)

        # correct the boundary condition numbers:
        # the index of the elements and neighbours have changed
        for iel, ibc, iface in product(
            range(nel1, self.nel), range(other.nbc), range(6)
        ):
            self.elem[iel].bcs[ibc, iface][1] = iel + 1
            bc = self.elem[iel].bcs[ibc, iface][0]
            if bc == "E" or bc == "P":
                neighbour = self.elem[iel].bcs[ibc, iface][3]
                self.elem[iel].bcs[ibc, iface][3] = neighbour + nel1

        # glue common faces together
        # only look for the neighbour in the first BC field because it should be the same in all fields.
        # FIXME: this will fail to correct periodic conditions if periodic domains are merged together.
        nfaces = 2 * self.ndim
        nchanges = 0  # counter for the boundary conditions connected
        if nbc == 0 or ignore_all_bcs:
            # Quickly exit the function
            logger.debug("no pairs of faces to merge")
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
                        dist2 = (x1 - x0) ** 2 + (y1 - y0) ** 2 + (z1 - z0) ** 2
                        if dist2 <= tol ** 2:
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
