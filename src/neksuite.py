#=============================================================================#
# neksuite                                                                    #
#                                                                             #
# A python module for reading and writing nek5000 files                       #
#                                                                             #
# Authors: Jacopo Canton, Nicolo' Fabbiane                                    #
# Contacts: jcanton(at)mech.kth.se, nicolo(at)mech.kth.se                     #
# Last edit: 2015-10-19                                                       #
#=============================================================================#
import struct
import numpy as np
from tvtk.api import tvtk


#==============================================================================
class datalims:
	"""
	    datalims
	    A class containing the extrema of all quantities stored in the mesh
	"""
	def __init__(self, var):
		#                    x,y,z   min,max
		self.pos  = np.zeros((3,      2))
		#self.pos  = np.zeros((var[0], 2))
		#                    u,v,w   min,max
		self.vel  = np.zeros((3,      2))
		#self.vel  = np.zeros((var[1], 2))
		#                    p       min,max
		self.pres = np.zeros((var[2], 2))
		#                    T       min,max
		self.temp = np.zeros((var[3], 2))
		#                    s_i     min,max
		self.scal = np.zeros((var[4], 2))

#==============================================================================
class nekelem:
	"""
	    nekelem
	    A class containing one nek element
	"""
	def __init__(self, var, lr1):
		#                    x,y,z   lz      ly      lx
		self.pos  = np.zeros((3,      lr1[2], lr1[1], lr1[0]))
		#                    u,v,w   lz      ly      lx
		self.vel  = np.zeros((3,      lr1[2], lr1[1], lr1[0]))
      #                    p       lz      ly      lx     
		self.pres = np.zeros((var[2], lr1[2], lr1[1], lr1[0]))
      #                    T       lz      ly      lx     
		self.temp = np.zeros((var[3], lr1[2], lr1[1], lr1[0]))
      #                    s_i     lz      ly      lx     
		self.scal = np.zeros((var[4], lr1[2], lr1[1], lr1[0]))

#==============================================================================
class nekdata:
	"""
	    nekdata
	    A class containing data for reading/writing binary nek files
	"""
	def __init__(self, ndim, nel, lr1, var):
		self.ndim   = ndim
		self.nel    = nel
		self.var    = var
		self.lr1    = lr1
		self.time   = []
		self.istep  = []
		self.wdsz   = []
		self.endian = []
		self.lims   = datalims(var)
		self.elem   = [nekelem(var, lr1) for i in range(nel)]


#==============================================================================
def readnek(fname):
	"""
	    readnek
	    A function for reading binary data from the nek5000 binary format

	    input variable:
	    fname : file name
	"""
	#
	try:
		infile = open(fname, 'rb')
	except IOError as e:
		print('I/O error ({0}): {1}'.format(e.errno, e.strerror))
		return -1
	#
	#---------------------------------------------------------------------------
	# READ HEADER
	#---------------------------------------------------------------------------
	#
	# read header
	header = infile.read(132).split()
	#
	# get word size
	wdsz = int(header[1])
	if (wdsz == 4):
		realtype = 'f'
	elif (wdsz == 8):
		realtype = 'd'
	else:
		print('ERROR: could not interpret real type (wdsz = %i)' %(wdsz))
		return -2
	#
	# get polynomial order
	lr1 = [int(header[2]), 
	       int(header[3]),
	       int(header[4])]
	#
	# compute total number of points per element
	npel = lr1[0] * lr1[1] * lr1[2]
	#
	# get number of pysical dimensions
	ndim = 2 + (lr1[2]>1)
	#
	# get number of elements
	nel = int(header[5])
	#
	# get number of elements in the file
	nelf = int(header[6])
	#
	# get current time
	time = float(header[7])
	#
	# get current time step
	istep = int(header[8])
	#
	# get file id
	fid = int(header[9])
	#
	# get tot number of files
	nf = int(header[10])
	#
	# get variables [XUPT]
	vars = header[-1]
	var = [0 for i in range(5)]
	for v in vars:
		if (v == 'X'):
			var[0] = ndim
		elif (v == 'U'):
			var[1] = ndim
		elif (v == 'P'):
			var[2] = 1
		elif (v == 'T'):
			var[3] = 1
		elif (v == 'S'):
			var[4] = 1 # TODO: need to know how this works
	#
	# compute number of scalar fields
	nfields = sum(var)
	#
	# identify endian encoding
	etagb = infile.read(4)
	etagL = struct.unpack('<f', etagb)[0]; etagL = int(etagL*1e5)/1e5
	etagB = struct.unpack('>f', etagb)[0]; etagB = int(etagB*1e5)/1e5
	if (etagL == 6.54321):
		# print('Reading little-endian file\n')
		emode = '<'
	elif (etagB == 6.54321):
		# print('Reading big-endian file\n')
		emode = '>'
	else:
		print('ERROR: could not interpret endianness')
		return -3
	#
	# read element map for the file
	elmap = infile.read(4*nelf)
	elmap = list(struct.unpack(emode+nelf*'i', elmap))
	#
	#---------------------------------------------------------------------------
	# READ DATA
	#---------------------------------------------------------------------------
	#
	# initialize data structure
	data = nekdata(ndim, nel, lr1, var)
	data.time   = time
	data.istep  = istep
	data.wdsz   = wdsz
	if (emode == '<'):
		data.endian = 'little'
	elif (emode == '>'):
		data.endian = 'big'
	#
	# read geometry
	data.lims.pos[:,0] =  float('inf')
	data.lims.pos[:,1] = -float('inf')
	for iel in elmap:
		for idim in range(var[0]): # if var[0] == 0, geometry is not read
			fi = infile.read(npel*wdsz)
			fi = list(struct.unpack(emode+npel*realtype, fi))
			ip = 0
			for iz in range(lr1[2]):
				for iy in range(lr1[1]):
					data.elem[iel-1].pos[idim,iz,iy,:] = fi[ip:ip+lr1[0]]
					ip += lr1[0]
			data.lims.pos[idim,0] = min([data.lims.pos[idim,0]]+fi)
			data.lims.pos[idim,1] = max([data.lims.pos[idim,1]]+fi)
	#
	# read velocity
	data.lims.vel[:,0] =  float('inf')
	data.lims.vel[:,1] = -float('inf')
	for iel in elmap:
		for idim in range(var[1]): # if var[1] == 0, velocity is not read
			fi = infile.read(npel*wdsz)
			fi = list(struct.unpack(emode+npel*realtype, fi))
			ip = 0
			for iz in range(lr1[2]):
				for iy in range(lr1[1]):
					data.elem[iel-1].vel[idim,iz,iy,:] = fi[ip:ip+lr1[0]]
					ip += lr1[0]
			data.lims.vel[idim,0] = min([data.lims.vel[idim,0]]+fi)
			data.lims.vel[idim,1] = max([data.lims.vel[idim,1]]+fi)
	#
	# read pressure 
	data.lims.pres[:,0] =  float('inf')
	data.lims.pres[:,1] = -float('inf')
	for iel in elmap:
		for ivar in range(var[2]): # if var[2] == 0, pressure is not read
			fi = infile.read(npel*wdsz)
			fi = list(struct.unpack(emode+npel*realtype, fi))
			ip = 0
			for iz in range(lr1[2]):
				for iy in range(lr1[1]):
					data.elem[iel-1].pres[ivar,iz,iy,:] = fi[ip:ip+lr1[0]]
					ip += lr1[0]
			data.lims.pres[ivar,0] = min([data.lims.pres[ivar,0]]+fi)
			data.lims.pres[ivar,1] = max([data.lims.pres[ivar,1]]+fi)
	#
	# read temperature
	data.lims.temp[:,0] =  float('inf')
	data.lims.temp[:,1] = -float('inf')
	for iel in elmap:
		for ivar in range(var[3]): # if var[3] == 0, temperature is not read
			fi = infile.read(npel*wdsz)
			fi = list(struct.unpack(emode+npel*realtype, fi))
			ip = 0
			for iz in range(lr1[2]):
				for iy in range(lr1[1]):
					data.elem[iel-1].temp[ivar,iz,iy,:] = fi[ip:ip+lr1[0]]
					ip += lr1[0]
			data.lims.temp[ivar,0] = min([data.lims.temp[ivar,0]]+fi)
			data.lims.temp[ivar,1] = max([data.lims.temp[ivar,1]]+fi)
	#
	# read scalar fields
	data.lims.scal[:,0] =  float('inf')
	data.lims.scal[:,1] = -float('inf')
	for iel in elmap:
		for ivar in range(var[4]): # if var[4] == 0, scalars are not read
			fi = infile.read(npel*wdsz)
			fi = list(struct.unpack(emode+npel*realtype, fi))
			ip = 0
			for iz in range(lr1[2]):
				for iy in range(lr1[1]):
					data.elem[iel-1].scal[ivar,iz,iy,:] = fi[ip:ip+lr1[0]]
					ip += lr1[0]
			data.lims.scal[ivar,0] = min([data.lims.scal[ivar,0]]+fi)
			data.lims.scal[ivar,1] = max([data.lims.scal[ivar,1]]+fi)
	#
	#
	# close file
	infile.close()
	#
	# output
	return data

#==============================================================================
def nek2vtk(field, downsample):
	"""
	    nek2vtk
	    A function for converting nek data to vtk data

	    input variables:
		 field : a dataset in nekdata format
		 downsample : flag T/F
	"""
	#
	if (downsample):
		ixs = field.lr1[0]-1
		iys = field.lr1[1]-1
		izs = max(field.lr1[2]-1,1)
	else:
		ixs = 1
		iys = 1
		izs = 1
	#
	iix = range(0,field.lr1[0],ixs); nix = len(iix)
	iiy = range(0,field.lr1[1],iys); niy = len(iiy)
	iiz = range(0,field.lr1[2],izs); niz = len(iiz)
	#
	nppel = nix*niy*niz
	nepel = (nix-1)*(niy-1)*max((niz-1),1)
	nel   = field.nel*nepel
	#
	if (field.ndim==3):
		nvert = 8
		cellType = tvtk.Hexahedron().cell_type
	else:
		nvert = 4
		cellType = tvtk.Quad().cell_type
	#
	ct = np.array(nel*[cellType])
	of = np.arange(0,nvert*nel,nvert)
	ce = np.zeros(nel*(nvert+1)); ce[range(0,nel*(nvert+1),nvert+1)] = nvert
	if (field.var[0]!=0):
		r  = np.zeros((nvert*nel,3))
	if (field.var[1]!=0):
		v  = np.zeros((nvert*nel,3))
	if (field.var[2]==1):
		p  = np.zeros((nvert*nel))
	if (field.var[3]==1):
		T  = np.zeros((nvert*nel))
	if (field.var[4]!=0):
		S  = np.zeros((nvert*nel,field.var[4]))
	#
	ice = -(nvert + 1)
	for iel in range(field.nel):
		for iz in range(niz):
			for iy in range(niy):
				for ix in range(nix):
					if (field.var[0]==3):
						r[iel*nppel + ix + iy*nix + iz*nix*niy,:] = field.elem[iel].pos [:, iiz[iz], iiy[iy], iix[ix]]
					if (field.var[1]==3):
						v[iel*nppel + ix + iy*nix + iz*nix*niy,:] = field.elem[iel].vel [:, iiz[iz], iiy[iy], iix[ix]]
					if (field.var[2]==1):
						p[iel*nppel + ix + iy*nix + iz*nix*niy]   = field.elem[iel].pres[:, iiz[iz], iiy[iy], iix[ix]]
					if (field.var[3]==1):
						T[iel*nppel + ix + iy*nix + iz*nix*niy]   = field.elem[iel].temp[:, iiz[iz], iiy[iy], iix[ix]]
					if (field.var[4]!=0):
						S[iel*nppel + ix + iy*nix + iz*nix*niy,:] = field.elem[iel].scal[:, iiz[iz], iiy[iy], iix[ix]]
		if (field.var[0]==3):
			for iz in max(range(niz-1), [0]):
				for iy in range(niy-1):
					for ix in range(nix-1):
						ice = ice + nvert + 1
						for face in range(field.ndim-1):
							ce[ice + face*4 + 1] = iel*nppel + ix   +  iy   *nix + (iz+face)*nix*niy
							ce[ice + face*4 + 2] = iel*nppel + ix+1 +  iy   *nix + (iz+face)*nix*niy
							ce[ice + face*4 + 3] = iel*nppel + ix+1 + (iy+1)*nix + (iz+face)*nix*niy
							ce[ice + face*4 + 4] = iel*nppel + ix   + (iy+1)*nix + (iz+face)*nix*niy
	
	# create the array of cells
	ca = tvtk.CellArray()
	ca.set_cells(nel, ce)
	# create the unstructured dataset
	dataset = tvtk.UnstructuredGrid(points=r)
	# set the cell types
	dataset.set_cells(ct, of, ca)
	# set the data
	dataset.point_data.vectors = v
	dataset.point_data.vectors.name = 'vel'
	if (field.var[2]==1):
		dataset.point_data.scalars = p
		dataset.point_data.scalars.name = 'pres'
	if (field.var[3]==1):
		dataset.point_data.add_array(T)
		dataset.point_data.get_array(2).name = 'temp'
	if (field.var[4]!=0):
		for ii in range(field.var[4]):
			dataset.point_data.add_array(S[:,ii])
			dataset.point_data.get_array(ii+3).name = 'scal_%d' % (ii+1)
	#
	dataset.point_data.update()
	#
	return dataset
