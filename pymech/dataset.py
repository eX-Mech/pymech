"""Interface for reading files as xarray datasets."""
from pathlib import Path
import numpy as np
import xarray as xr

from .neksuite import readnek


def open_dataset(path, chunks=None, parallel=False):
	"""Helper function for opening a file as an xarray dataset."""
	path = Path(path)
	if path.suffix.startswith('.f'):
		_open = _open_nek_dataset
	else:
		raise NotImplementedError(
			"Filetype: {} is not supported.".format(path.suffix)
		)

	return _open(path, chunks, parallel)


#  @profile
def _open_nek_dataset(path, chunks, parallel):
	"""Interface for converting Nek field files into xarray datasets

	If parallel
	"""

	field = readnek(path)
	if isinstance(field, int):
		raise IOError("Failed to load {}".format(path))

	elements = field.elem

	if parallel:
		import dask
		import dask.bag as db

		#
		elements = db.from_sequence(elements, npartitions=4)

		DataStore = dask.delayed(_NekDataStore)
		create_dset = dask.delayed(_create_nek_dataset)
		set_coords = dask.delayed(_set_nek_coords)

		def tasks(elem):
			return set_coords(create_dset(DataStore(elem)))

		task_graph = elements.map(tasks)

		#  dask.visualize(task_graph, filename="graph.svg")
		elem_dsets = dask.compute(*task_graph)
	else:
		DataStore = _NekDataStore
		create_dset = _create_nek_dataset
		set_coords = _set_nek_coords

		elem_dsets = [
			set_coords(create_dset(DataStore(elem), chunks))
			for elem in elements
		]


	# See: https://github.com/MITgcm/xmitgcm/pull/200
	if xr.__version__ < '0.15.2':
		ds = xr.combine_by_coords(elem_dsets)
	else:
		ds = xr.combine_by_coords(elem_dsets, combine_attrs="drop")

	ds.coords.update({"time": field.time})

	return ds


def _create_nek_dataset(store, chunks=False):
	ds = xr.Dataset.load_store(store)

	if chunks:
		from dask.base import tokenize

		elem_id = id(store)
		token = tokenize(chunks, elem_id)
		name_prefix = "create_nek_dataset-%d" % token

		return ds.chunk(chunks, name_prefix=name_prefix, token=token)
	else:
		return ds


def _set_nek_coords(ds):
	axes = ("z", "y", "x")
	return ds.set_coords(axes)


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
