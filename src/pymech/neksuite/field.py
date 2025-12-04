from __future__ import annotations

import io
import os
import struct
import sys
from pathlib import Path
from typing import BinaryIO, Optional, Tuple, Union

import numpy as np
from attrs import define, field

from pymech.core import HexaData
from pymech.log import logger


def _as_unicode(string: StringOrBytes) -> str:
    if isinstance(string, bytes):
        return string.decode()
    else:
        return string


def _as_tuple_of_ints(seq):
    return tuple(int(s) for s in seq)


StringOrBytes = Union[str, bytes]
PathLike = Union[str, os.PathLike]


@define
class Header:
    """Dataclass for Nek5000 field file header. This relies on the package
    attrs_ and its ability to do type-validation and type-conversion of the
    header metadata.

    .. _attrs: https://www.attrs.org/en/stable/

    """

    # get word size: single or double precision
    wdsz: int = field(converter=int)
    # get polynomial order
    orders: Tuple[int, ...] = field(converter=_as_tuple_of_ints)
    # get number of elements
    nb_elems: int = field(converter=int)
    # get number of elements in the file
    nb_elems_file: int = field(converter=int)
    # get current time
    time: float = field(converter=float)
    # get current time step
    istep: int = field(converter=int)
    # get file id
    fid: int = field(converter=int)
    # get tot number of files
    nb_files: int = field(converter=int)

    # NOTE: field(factory=...) specifies the default value for the field.
    # https://www.attrs.org/en/stable/init.html#defaults

    # get variables [XUPTS[01-99]]
    variables: str = field(converter=_as_unicode, factory=str)
    # floating point precision
    realtype: str = field(factory=str)
    # compute total number of points per element
    nb_pts_elem: int = field(factory=int)
    # get number of physical dimensions
    nb_dims: int = field(factory=int)
    # get number of variables
    nb_vars: Tuple[int, ...] = field(factory=tuple)

    def __attrs_post_init__(self):
        # get word size: single or double precision
        wdsz = self.wdsz
        if not self.realtype:
            if wdsz == 4:
                self.realtype = "f"
            elif wdsz == 8:
                self.realtype = "d"
            else:
                raise ValueError(f"Could not interpret real type (wdsz = {wdsz})")

        orders = self.orders
        if not self.nb_pts_elem:
            self.nb_pts_elem = np.prod(orders)

        if not self.nb_dims:
            self.nb_dims = 2 + int(orders[2] > 1)

        if not self.variables and not self.nb_vars:
            raise ValueError("Both variables and nb_vars cannot be uninitialized.")
        elif self.variables:
            self.nb_vars = self._variables_to_nb_vars()
        elif self.nb_vars:
            self.variables = self._nb_vars_to_variables()

        logger.debug(f"Variables: {self.variables}, nb_vars: {self.nb_vars}")

    def _variables_to_nb_vars(self) -> Optional[Tuple[int, ...]]:
        # get variables [XUPTS[01-99]]
        variables = self.variables
        nb_dims = self.nb_dims

        if not variables:
            raise ValueError("Failed to convert variables to nb_vars")

        if not nb_dims:
            raise ValueError("Unintialized nb_dims")

        def nb_scalars():
            index_s = variables.index("S")
            return int(variables[index_s + 1 :])

        nb_vars = (
            nb_dims if "X" in variables else 0,
            nb_dims if "U" in variables else 0,
            1 if "P" in variables else 0,
            1 if "T" in variables else 0,
            nb_scalars() if "S" in variables else 0,
        )

        return nb_vars

    def _nb_vars_to_variables(self) -> Optional[str]:
        nb_vars = self.nb_vars
        if not nb_vars:
            raise ValueError("Failed to convert nb_vars to variables")

        str_vars = ("X", "U", "P", "T", f"S{nb_vars[4]:02d}")
        variables = (str_vars[i] if nb_vars[i] > 0 else "" for i in range(5))
        return "".join(variables)

    def as_bytestring(self) -> bytes:
        header = "#std %1i %2i %2i %2i %10i %10i %20.13E %9i %6i %6i %s" % (
            self.wdsz,
            self.orders[0],
            self.orders[1],
            self.orders[2],
            self.nb_elems,
            self.nb_elems_file,
            self.time,
            self.istep,
            self.fid,
            self.nb_files,
            self.variables,
        )
        return header.ljust(132).encode("utf-8")


def read_header(path_or_file_obj: Union[PathLike, BinaryIO]) -> Header:
    """Make a :class:`pymech.neksuite.Header` instance from a file buffer
    opened in binary mode.

    """
    if isinstance(path_or_file_obj, (str, os.PathLike)):
        with Path(path_or_file_obj).open("rb") as fp:
            header = fp.read(132).split()
    elif isinstance(path_or_file_obj, io.BufferedReader):
        fp = path_or_file_obj
        header = fp.read(132).split()
    else:
        raise ValueError("Should be a path or opened file object in 'rb' mode.")

    logger.debug(b"Header: " + b" ".join(header))
    if len(header) < 12:
        raise IOError("Header of the file was too short.")

    # Relying on attrs converter to type-cast. Mypy will complain
    return Header(header[1], header[2:5], *header[5:12])  # type: ignore[arg-type]


# ==============================================================================
def readnek(fname, dtype="float64", skip_vars=()):
    """A function for reading binary data from the nek5000 binary format

    Parameters
    ----------
    fname : str
        File name
    dtype : str or type
        Floating point data type. See also :class:`pymech.core.Elem`.
    skip_vars: tuple[str]
        Variables to skip. Valid values to skip are ``("x", "y", "z", "ux",
        "uy", "uz", "pressure", "temperature", "s01", "s02", ...)``.  It also
        accept some extra values ``("vx", "vy", "vz", "p", "t")``.  If empty
        (default), it reads all variables available in the file.

    """
    #
    infile = open(fname, "rb")
    #
    # ---------------------------------------------------------------------------
    # READ HEADER
    # ---------------------------------------------------------------------------
    #
    # read header
    h = read_header(infile)
    #
    # identify endian encoding
    etagb = infile.read(4)
    etagL = struct.unpack("<f", etagb)[0]
    etagL = int(etagL * 1e5) / 1e5
    etagB = struct.unpack(">f", etagb)[0]
    etagB = int(etagB * 1e5) / 1e5
    if etagL == 6.54321:
        logger.debug("Reading little-endian file\n")
        emode = "<"
    elif etagB == 6.54321:
        logger.debug("Reading big-endian file\n")
        emode = ">"
    else:
        raise ValueError("Could not interpret endianness")

    #
    # read element map for the file
    elmap = infile.read(4 * h.nb_elems_file)
    elmap = struct.unpack(emode + h.nb_elems_file * "i", elmap)
    #
    # ---------------------------------------------------------------------------
    # READ DATA
    # ---------------------------------------------------------------------------
    #
    # initialize data structure
    data = HexaData(h.nb_dims, h.nb_elems, h.orders, h.nb_vars, 0, dtype)
    data.time = h.time
    data.istep = h.istep
    data.wdsz = h.wdsz
    data.elmap = np.array(elmap, dtype=np.int32)
    if emode == "<":
        data.endian = "little"
    elif emode == ">":
        data.endian = "big"

    bytes_elem = h.nb_pts_elem * h.wdsz

    def read_file_into_data(data_var, index_var):
        """Read binary file into an array attribute of ``data.elem``"""
        fi = infile.read(bytes_elem)
        fi = np.frombuffer(fi, dtype=emode + h.realtype, count=h.nb_pts_elem)

        # Replace elem array in-place with
        # array read from file after reshaping as
        elem_shape = h.orders[::-1]  # lz, ly, lx
        data_var[index_var, ...] = fi.reshape(elem_shape)

    def skip_elements(nb_elements=1):
        infile.seek(bytes_elem * nb_elements, os.SEEK_CUR)

    # read geometry
    geometry_vars = "x", "y", "z"
    nb_vars = h.nb_vars[0]
    skip_condition = tuple(geometry_vars[idim] in skip_vars for idim in range(nb_vars))
    if nb_vars:
        if all(skip_condition):
            skip_elements(h.nb_elems * nb_vars)
        else:
            if 0 in elmap:
                logger.warning(
                    "The 'elmap' appears to be corrupted as it contains an unexpected zero value."
                    " As a workaround, Pymech will read data by iterating over the entire set of"
                    " elements instead of following the map provided by 'elmap'."
                )
                element_idxs = range(h.nb_elems_file)
            else:
                element_idxs = (idx - 1 for idx in elmap)

            for iel in element_idxs:
                el = data.elem[iel]
                for idim in range(nb_vars):
                    if skip_condition[idim]:
                        skip_elements()
                    else:
                        read_file_into_data(el.pos, idim)

    # read velocity
    velocity_vars1 = "ux", "uy", "uz"
    velocity_vars2 = "vx", "vy", "vz"
    nb_vars = h.nb_vars[1]
    skip_condition1 = tuple(
        velocity_vars1[idim] in skip_vars for idim in range(nb_vars)
    )
    skip_condition2 = tuple(
        velocity_vars2[idim] in skip_vars for idim in range(nb_vars)
    )

    if nb_vars:
        if all(skip_condition1) or all(skip_condition2):
            skip_elements(h.nb_elems * nb_vars)
        else:
            if 0 in elmap:
                element_idxs = range(h.nb_elems_file)
            else:
                element_idxs = (idx - 1 for idx in elmap)

            for iel in element_idxs:
                el = data.elem[iel]
                for idim in range(nb_vars):
                    if skip_condition1[idim] or skip_condition2[idim]:
                        skip_elements()
                    else:
                        read_file_into_data(el.vel, idim)

    #
    # read pressure
    nb_vars = h.nb_vars[2]
    skip_condition = any({"p", "pressure"}.intersection(skip_vars))
    if nb_vars:
        if skip_condition:
            skip_elements(h.nb_elems * nb_vars)
        else:
            if 0 in elmap:
                element_idxs = range(h.nb_elems_file)
            else:
                element_idxs = (idx - 1 for idx in elmap)

            for iel in element_idxs:
                el = data.elem[iel]
                for ivar in range(nb_vars):
                    read_file_into_data(el.pres, ivar)

    #
    # read temperature
    nb_vars = h.nb_vars[3]
    skip_condition = any({"t", "temperature"}.intersection(skip_vars))
    if nb_vars:
        if skip_condition:
            skip_elements(h.nb_elems * nb_vars)
        else:
            if 0 in elmap:
                element_idxs = range(h.nb_elems_file)
            else:
                element_idxs = (idx - 1 for idx in elmap)

            for iel in element_idxs:
                el = data.elem[iel]
                for ivar in range(nb_vars):
                    read_file_into_data(el.temp, ivar)

    #
    # read scalar fields
    #
    nb_vars = h.nb_vars[4]
    scalar_vars = tuple(f"s{i:02d}" for i in range(1, nb_vars + 1))
    skip_condition = tuple(scalar_vars[ivar] in skip_vars for ivar in range(nb_vars))
    if nb_vars:
        if all(skip_condition):
            skip_elements(h.nb_elems * nb_vars)
        else:
            # NOTE: This is not a bug!
            # Unlike other variables, scalars are in the outer loop and elements
            # are in the inner loop
            for ivar in range(nb_vars):
                if skip_condition[ivar]:
                    skip_elements(h.nb_elems)
                else:
                    if 0 in elmap:
                        element_idxs = range(h.nb_elems_file)
                    else:
                        element_idxs = (idx - 1 for idx in elmap)

                    for iel in element_idxs:
                        el = data.elem[iel]
                        read_file_into_data(el.scal, ivar)

    #
    #
    # close file
    infile.close()
    #
    # output
    return data


# ==============================================================================
def writenek(fname, data):
    """A function for writing binary data in the nek5000 binary format

    Parameters
    ----------
    fname : str
            file name
    data : :class:`pymech.core.HexaData`
            data structure
    """
    #
    outfile = open(fname, "wb")
    #
    # ---------------------------------------------------------------------------
    # WRITE HEADER
    # ---------------------------------------------------------------------------
    #
    h = Header(
        data.wdsz,
        data.lr1,
        data.nel,
        data.nel,
        data.time,
        data.istep,
        fid=0,
        nb_files=1,
        nb_vars=data.var,
    )
    # NOTE: multiple files (not implemented). See fid, nb_files, nb_elem_file above
    #
    # get fields to be written
    #
    # get word size
    if h.wdsz == 4:
        logger.debug("Writing single-precision file")
    elif h.wdsz == 8:
        logger.debug("Writing double-precision file")
    else:
        raise ValueError("Could not interpret real type (wdsz = %i)" % (data.wdsz))
    #
    # generate header
    outfile.write(h.as_bytestring())
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
    # write element map for the file
    correct_endianness(data.elmap).tofile(outfile)
    #
    # ---------------------------------------------------------------------------
    # WRITE DATA
    # ---------------------------------------------------------------------------
    #
    # compute total number of points per element
    #  npel = data.lr1[0] * data.lr1[1] * data.lr1[2]

    def write_ndarray_to_file(a):
        """Write a data array to the output file in the requested precision and endianness"""
        if data.wdsz == 4:
            correct_endianness(a.astype(np.float32)).tofile(outfile)
        else:
            correct_endianness(a).tofile(outfile)

    #
    # write geometry
    for iel in data.elmap:
        for idim in range(data.var[0]):  # if var[0] == 0, geometry is not written
            write_ndarray_to_file(data.elem[iel - 1].pos[idim, :, :, :])
    #
    # write velocity
    for iel in data.elmap:
        for idim in range(data.var[1]):  # if var[1] == 0, velocity is not written
            write_ndarray_to_file(data.elem[iel - 1].vel[idim, :, :, :])
    #
    # write pressure
    for iel in data.elmap:
        for ivar in range(data.var[2]):  # if var[2] == 0, pressure is not written
            write_ndarray_to_file(data.elem[iel - 1].pres[ivar, :, :, :])
    #
    # write temperature
    for iel in data.elmap:
        for ivar in range(data.var[3]):  # if var[3] == 0, temperature is not written
            write_ndarray_to_file(data.elem[iel - 1].temp[ivar, :, :, :])
    #
    # write scalars
    #
    # NOTE: This is not a bug!
    # Unlike other variables, scalars are in the outer loop and elements
    # are in the inner loop
    #
    for ivar in range(data.var[4]):  # if var[4] == 0, scalars are not written
        for iel in data.elmap:
            write_ndarray_to_file(data.elem[iel - 1].scal[ivar, :, :, :])
    #
    # write max and min of every field in every element (forced to single precision)
    if data.ndim == 3:
        #
        for iel in data.elmap:
            for idim in range(data.var[0]):
                correct_endianness(
                    np.min(data.elem[iel - 1].pos[idim, :, :, :]).astype(np.float32)
                ).tofile(outfile)
                correct_endianness(
                    np.max(data.elem[iel - 1].pos[idim, :, :, :]).astype(np.float32)
                ).tofile(outfile)
        for iel in data.elmap:
            for idim in range(data.var[1]):
                correct_endianness(
                    np.min(data.elem[iel - 1].vel[idim, :, :, :]).astype(np.float32)
                ).tofile(outfile)
                correct_endianness(
                    np.max(data.elem[iel - 1].vel[idim, :, :, :]).astype(np.float32)
                ).tofile(outfile)
        for iel in data.elmap:
            for ivar in range(data.var[2]):
                correct_endianness(
                    np.min(data.elem[iel - 1].pres[ivar, :, :, :]).astype(np.float32)
                ).tofile(outfile)
                correct_endianness(
                    np.max(data.elem[iel - 1].pres[ivar, :, :, :]).astype(np.float32)
                ).tofile(outfile)
        for iel in data.elmap:
            for ivar in range(data.var[3]):
                correct_endianness(
                    np.min(data.elem[iel - 1].temp[ivar, :, :, :]).astype(np.float32)
                ).tofile(outfile)
                correct_endianness(
                    np.max(data.elem[iel - 1].temp[ivar, :, :, :]).astype(np.float32)
                ).tofile(outfile)
        for iel in data.elmap:
            for ivar in range(data.var[4]):
                correct_endianness(
                    np.min(data.elem[iel - 1].scal[ivar, :, :, :]).astype(np.float32)
                ).tofile(outfile)
                correct_endianness(
                    np.max(data.elem[iel - 1].scal[ivar, :, :, :]).astype(np.float32)
                ).tofile(outfile)

    # close file
    outfile.close()
    #
    # output
    return 0
