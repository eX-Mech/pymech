"""Interface for reading files as xarray datasets."""
from pathlib import Path
import numpy as np
import xarray as xr

from .neksuite import readnek


def open_dataset(path):
	"""Helper function for opening a file as an xarray dataset."""
	path = Path(path)
	if path.suffix.startswith('.f'):
		_open = open_nek_dataset
	else:
		raise NotImplementedError(
			"Filetype: {} is not supported.".format(path.suffix)
		)

	return _open(path)


#  @profile
def open_nek_dataset(path):
	"""Interface for converting Nek field files into xarray datasets"""
	field = readnek(path)
	if isinstance(field, int):
		raise IOError("Failed to load {}".format(path))

	elements = field.elem
	elem_stores = [_NekDataStore(elem) for elem in elements]
	elem_dsets = [
		xr.Dataset.load_store(store).set_coords(store.axes)
		for store in elem_stores
	]

	# See: https://github.com/MITgcm/xmitgcm/pull/200
	if xr.__version__ < '0.15.2':
		ds = xr.combine_by_coords(elem_dsets)
	else:
		ds = xr.combine_by_coords(elem_dsets, combine_attrs="drop")

	ds.coords.update({"time": field.time})

	return ds


class _NekDataStore(xr.backends.common.AbstractDataStore):
	"""Xarray store for a Nek field element."""
	def __init__(self, elem):
		self.elem = elem
		self.axes = ("z", "y", "x")

	def meshgrid_to_dim(self, mesh):
		"""Reverse of np.meshgrid. This method extracts one-dimensional
		coordinates from a cubical array format for every direction

		"""
		dim = np.unique(np.round(mesh, 8))
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
		return attrs

	#  @profile
	def get_variables(self):
		"""Generate an xarray dataset from a single element."""
		ax = self.axes
		elem = self.elem

		data_vars = {
			ax[2]: self.meshgrid_to_dim(elem.pos[0]),  # x
			ax[1]: self.meshgrid_to_dim(elem.pos[1]),  # y
			ax[0]: self.meshgrid_to_dim(elem.pos[2]),  # z
			"xmesh": (ax, elem.pos[0]),
			"ymesh": (ax, elem.pos[1]),
			"zmesh": (ax, elem.pos[2]),
			"ux": (ax, elem.vel[0]),
			"uy": (ax, elem.vel[1]),
			"uz": (ax, elem.vel[2]),
		}
		if elem.pres.size:
			data_vars["pressure"] = ax, elem.pres[0]

		if elem.temp.size:
			data_vars["temperature"] = ax, elem.temp[0]

		if elem.scal.size:
			data_vars.update(
				{
					"s{:02d}".format(iscalar + 1): (ax, elem.scal[iscalar])
					for iscalar in range(elem.scal.shape[0])
				}
			)

		return data_vars
