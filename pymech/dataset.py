"""Interface for reading files as xarray datasets."""
from pathlib import Path
import numpy as np
import xarray as xr

from .neksuite import readnek


def open_dataset(path):
	"""Helper function for opening a file as an xarray dataset."""
	path = Path(path)
	if path.suffix.startswith('.f'):
		cls = NekDataset
	else:
		raise NotImplementedError(
			"Filetype: {} is not supported.".format(path.suffix)
		)

	return cls(path).open()


class NekDataset:
	"""Interface for converting Nek field files into xarray datasets"""
	def __init__(self, path):
		self.field = readnek(path)
		self.axes = ("z", "y", "x")

	def meshgrid_to_dim(self, mesh):
		"""Reverse of np.meshgrid. This method extracts one-dimensional
		coordinates from a cubical array format for every direction

		"""
		dim = np.unique(np.round(mesh, 8))
		return dim

	def elem_to_dataset(self, elem):
		"""Generate an xarray dataset from a single element."""
		ax = self.axes

		ds = xr.Dataset(
			data_vars={
				"ux": (ax, elem.vel[0]),
				"uy": (ax, elem.vel[1]),
				"uz": (ax, elem.vel[2]),
			},
			coords={
				ax[2]: self.meshgrid_to_dim(elem.pos[0]),  # x
				ax[1]: self.meshgrid_to_dim(elem.pos[1]),  # y
				ax[0]: self.meshgrid_to_dim(elem.pos[2]),  # z
				"xmesh": (ax, elem.pos[0]),
				"ymesh": (ax, elem.pos[1]),
				"zmesh": (ax, elem.pos[2]),
			},
			attrs={
				"boundary_conditions": elem.bcs,
				"curvature": elem.curv,
				"curvature_type": elem.ccurv,
			}
		)

		if elem.pres.size:
			ds["pressure"] = ax, elem.pres[0]

		if elem.temp.size:
			ds["temperature"] = ax, elem.temp[0]

		if elem.scal.size:
			ds["scalar"] = ax, elem.scal[0]

		return ds

	def open(self):
		elements = self.field.elem
		ds = xr.combine_by_coords(map(self.elem_to_dataset, elements))
		ds.coords.update({"time": self.field.time})
		return ds
