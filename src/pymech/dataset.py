import re
from functools import partial
from pathlib import Path

import numpy as np
import uxarray as uxr
import xarray as xr
from xarray.core.utils import Frozen

from .neksuite import readnek

__all__ = (
    "open_dataset",
    "open_mfdataset",
)


nek_ext_pattern = re.compile(
    r"""
   .*          # one or more characters
   \.          # character "."
   f           # character "f"
   (\d{5}|ld)  # 5 digits or the characters "ld"
""",
    re.VERBOSE,
)


def can_open_nek_dataset(path):
    """A regular expression check of the file extension.

    .. hint::

        - Would not match: .f90 .f .fort .f0000
        - Would match: .fld .f00001 .f12345

    """
    return nek_ext_pattern.match(str(path))


def open_dataset(path, **kwargs):
    """Helper function for opening a file as an :class:`xarray.Dataset`.

    Parameters
    ----------
    path : str
            Path to a field file (only Nek files are supported at the moment.)

    kwargs : dict
            Keyword arguments passed on to the compatible open function.

    """
    if can_open_nek_dataset(path):
        _open = _open_nek_dataset
    else:
        raise NotImplementedError(f"Filetype: {Path(path).suffix} is not supported.")

    return _open(path, **kwargs)


open_mfdataset = partial(
    xr.open_mfdataset, combine="nested", concat_dim="time", engine="pymech"
)
open_mfdataset.__doc__ = """Helper function for opening multiple files as an
:class:`xarray.Dataset`. See :func:`xarray.open_mfdataset` for documentation on
parameters."""


def _open_nek_dataset(path, drop_variables=None):
    """Interface for converting Nek field files into xarray_ datasets.

    .. _xarray: https://docs.xarray.dev/en/stable/
    """
    field = readnek(path)
    if isinstance(field, int):
        raise OSError(f"Failed to load {path}")

    elements = field.elem
    elem_stores = [_NekDataStore(elem) for elem in elements]
    try:
        elem_dsets = [
            xr.Dataset.load_store(store).set_coords(store.axes) for store in elem_stores
        ]
    except ValueError as err:
        raise NotImplementedError(
            "Opening dataset failed because you probably tried to open a field file "
            "with an unsupported mesh. "
            "The `pymech.open_dataset` function currently works only with cartesian "
            "box meshes. For more details on this, see "
            "https://github.com/eX-Mech/pymech/issues/31"
        ) from err

    # See: https://github.com/MITgcm/xmitgcm/pull/200
    ds = xr.combine_by_coords(elem_dsets, combine_attrs="drop")
    ds.coords.update({"time": field.time})

    if drop_variables:
        ds = ds.drop_vars(drop_variables)

    return ds


def extract_elem_data(elem_array):
    # Use lists to accumulate data
    x_list, y_list, z_list = [], [], []
    ux_list, uy_list, uz_list, p_list = [], [], [], []

    # Loop through the elements in the array
    for elem in elem_array:
        # Append data to respective lists
        x_list.append(np.ravel(elem.pos[0]))
        y_list.append(np.ravel(elem.pos[1]))
        z_list.append(np.ravel(elem.pos[2]))

        ux_list.append(np.ravel(elem.vel[0]))
        uy_list.append(np.ravel(elem.vel[1]))
        uz_list.append(np.ravel(elem.vel[2]))

        p_list.append(np.ravel(elem.pres))

    # Convert lists to NumPy arrays after the loop
    x = np.concatenate(x_list)
    y = np.concatenate(y_list)
    z = np.concatenate(z_list)
    ux = np.concatenate(ux_list)
    uy = np.concatenate(uy_list)
    uz = np.concatenate(uz_list)
    p = np.concatenate(p_list)

    # Create an xarray dataset
    data = xr.Dataset(
        {
            "ux": (["points"], ux),
            "uy": (["points"], uy),
            "uz": (["points"], uz),  # Correctly assign uz
            "p": (["points"], p),  # Correctly assign p
        },
        coords={
            "x": (["points"], x),
            "y": (["points"], y),
            "z": (["points"], z),
        },
    )

    # Wrap the xarray dataset in a uxarray grid (optional)
    ux_ds = uxr.Grid.from_dataset(data)

    return ux_ds


class PymechXarrayBackend(xr.backends.BackendEntrypoint):
    def guess_can_open(self, filename_or_obj):
        return can_open_nek_dataset(filename_or_obj)

    def open_dataset(
        self,
        filename_or_obj,
        *,
        drop_variables=None,
        # other backend specific keyword arguments
        # `chunks` and `cache` DO NOT go here, they are handled by xarray
    ):
        return _open_nek_dataset(filename_or_obj, drop_variables)

    open_dataset_parameters = ("filename_or_obj", "drop_variables")


class _NekDataStore(xr.backends.common.AbstractDataStore):
    """Xarray store for a Nek field element.

    Parameters
    ----------
    elem: :class:`pymech.core.Elem`
        A Nek5000 element.

    """

    axes = ("z", "y", "x")

    def __init__(self, elem):
        self.elem = elem

    def meshgrid_to_dim(self, mesh):
        """Reverse of np.meshgrid. This method extracts one-dimensional
        coordinates from a cubical array format for every direction

        """
        u, indices = np.unique(np.round(mesh, 8), return_index=True)
        mesh_1d = np.reshape(mesh, np.size(mesh))
        dim = mesh_1d[indices]
        return dim

    def get_dimensions(self):
        return self.axes

    def get_attrs(self):
        elem = self.elem
        attrs = {
            "boundary_conditions": elem.bcs,
            "curvature": elem.curv,
            "curvature_type": elem.ccurv,
        }
        return Frozen(attrs)

    def get_variables(self):
        """Generate an xarray dataset from a single element."""
        ax = self.axes
        elem = self.elem

        data_vars = {
            ax[2]: self.meshgrid_to_dim(elem.pos[0]),  # x
            ax[1]: self.meshgrid_to_dim(elem.pos[1]),  # y
            ax[0]: self.meshgrid_to_dim(elem.pos[2]),  # z
            "xmesh": xr.Variable(ax, elem.pos[0]),
            "ymesh": xr.Variable(ax, elem.pos[1]),
            "zmesh": xr.Variable(ax, elem.pos[2]),
            "ux": xr.Variable(ax, elem.vel[0]),
            "uy": xr.Variable(ax, elem.vel[1]),
            "uz": xr.Variable(ax, elem.vel[2]),
        }
        if elem.pres.size:
            data_vars["pressure"] = xr.Variable(ax, elem.pres[0])

        if elem.temp.size:
            data_vars["temperature"] = xr.Variable(ax, elem.temp[0])

        if elem.scal.size:
            data_vars.update(
                {
                    "s{:02d}".format(iscalar + 1): xr.Variable(ax, elem.scal[iscalar])
                    for iscalar in range(elem.scal.shape[0])
                }
            )

        return Frozen(data_vars)


def open_unstruc_dataset(path):

    # Proposed Methodology
    # Step 1: Use readnek to import the data
    # Step 2: Create an array of nodes, elements, and fields
    # Step 3: Create a grid and xarray dataset from the data array
    # Step 4: Create a uxarray dataset and return it

    field = readnek(path)
    if isinstance(field, int):
        raise OSError(f"Failed to load {path}")

    elements = field.elem

    # Method 1 : adapt the existing method used for xarray

    # elem_stores = [_NekDataStore(elem) for elem in elements]
    #
    # try:
    #     elem_dsets = [
    #         uxr.UxDataset.load_store(store).set_coords(store.axes) for store in elem_stores
    #     ]
    # except Exception as error:
    #     print("uxarray failure")
    #     print(error)
    #     print(elem_stores[0].axes)

    # Method 2 : manually create array of x, y, z and variables and use it to
    # make a uxarray dataset
    ds = extract_elem_data(elements)
