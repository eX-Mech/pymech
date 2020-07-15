"""Interface for reading files as xarray_ datasets.

.. _xarray: https://xarray.pydata.org

"""
from pathlib import Path
from glob import glob

import numpy as np
import xarray as xr

from .neksuite import readnek


def open_dataset(path, **kwargs):
	"""Helper function for opening a file as an xarray_ dataset.

	Parameters
	----------
	path : str
		Path to a field file (only Nek files are supported at the moment.)

	kwargs : dict
		Keyword arguments passed on to the compatible open function.

	"""
	path = Path(path)
	if path.suffix.startswith('.f'):
		_open = _open_nek_dataset
	else:
		raise NotImplementedError(
			"Filetype: {} is not supported.".format(path.suffix)
		)

	return _open(path, **kwargs)


def open_mfdataset(
	paths,
	chunks=None,
	concat_dim="time",
	compat="no_conflicts",
	preprocess=None,
	engine=None,
	lock=None,
	data_vars="all",
	coords="different",
	combine="nested",
	autoclose=None,
	parallel=False,
	join="outer",
	attrs_file=None,
	**kwargs,
):
	"""Helper function for opening multiple files as an xarray_ dataset.
	Adapted from upstream implementation_. See docs_ for usage.

	.. todo::

		To be removed when a backend entrypoint_ is implementated.

	.. _implementation: https://github.com/pydata/xarray/blob/484d1ce5ff8969b6ca6fa942b344379725f33b9c/xarray/backends/api.py#L726
	.. _docs: https://xarray.pydata.org/en/v0.15.1/generated/xarray.open_mfdataset.html
	.. _entrypoint: https://github.com/pydata/xarray/pull/3166

	"""
	if isinstance(paths, str):
		paths = sorted(glob(paths))
	else:
		paths = [str(p) if isinstance(p, Path) else p for p in paths]

	if not paths:
		raise OSError("no files to open")

	# If combine='by_coords' then this is unnecessary, but quick.
	# If combine='nested' then this creates a flat list which is easier to
	# iterate over, while saving the originally-supplied structure as "ids"
	if combine == "nested":
		if isinstance(concat_dim, (str, xr.DataArray)) or concat_dim is None:
			concat_dim = [concat_dim]

	open_kwargs = dict()

	if parallel:
		import dask

		# wrap the open_dataset, getattr, and preprocess with delayed
		open_ = dask.delayed(open_dataset)
		if preprocess is not None:
			preprocess = dask.delayed(preprocess)
	else:
		open_ = open_dataset

	datasets = [open_(p, **open_kwargs) for p in paths]
	if preprocess is not None:
		datasets = [preprocess(ds) for ds in datasets]

	if parallel:
		# calling compute here will return the datasets
		# the underlying datasets will still be stored as dask arrays
		datasets, = dask.compute(datasets)

	# Combine all datasets, closing them in case of a ValueError
	try:
		if combine == "nested":
			# Combined nested list by successive concat and merge operations
			# along each dimension, using structure given by "ids"
			combined = xr.combine_nested(
				datasets,
				concat_dim=concat_dim,
				compat=compat,
				data_vars=data_vars,
				coords=coords,
				join=join,
			)
		elif combine == "by_coords":
			# Redo ordering from coordinates, ignoring how they were ordered
			# previously
			combined = xr.combine_by_coords(
				datasets, compat=compat, data_vars=data_vars, coords=coords, join=join
			)
		else:
			raise ValueError(
				"{} is an invalid option for the keyword argument"
				" ``combine``".format(combine)
			)
	except ValueError:
		for ds in datasets:
			ds.close()
		raise

	# read global attributes from the attrs_file or from the first dataset
	if attrs_file is not None:
		if isinstance(attrs_file, Path):
			attrs_file = str(attrs_file)
		combined.attrs = datasets[paths.index(attrs_file)].attrs
	else:
		combined.attrs = datasets[0].attrs

	return combined


def _open_nek_dataset(path):
	"""Interface for converting Nek field files into xarray_ datasets."""
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
