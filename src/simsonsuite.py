#=============================================================================#
# simsonsuite                                                                 #
#                                                                             #
# A python module for reading and writing SIMSON files                        #
#                                                                             #
# Authors: Jacopo Canton, Nicolo' Fabbiane                                    #
# Contacts: jcanton(at)mech.kth.se, nicolo(at)mech.kth.se                     #
# Last edit: 2015-10-19                                                       #
#=============================================================================#
import struct
import numpy as np


#==============================================================================
class datalims:
	"""
	    datalims
	    A class containing the extrema of all quantities stored in the mesh
	"""
	def __init__(self, var):
		#                    x,y,z   min,max
		self.pos  = np.zeros((3     , 2))
		#                    u,v,w   min,max
		self.vel  = np.zeros((3     , 2))
		#                    p       min,max
		self.pres = np.zeros((var[2], 2))
		#                    T       min,max
		self.temp = np.zeros((var[3], 2))
		#                    s_i     min,max
		self.scal = np.zeros((var[4], 2))

#==============================================================================
class elem:
	"""
	    elem
	    A class containing one nek element/SIMSON flow field
	"""
	def __init__(self, var, lr1):
		#                    x,y,z   lz      ly      lx
		self.pos  = np.zeros((3     , lr1[2], lr1[1], lr1[0]))
		#                    u,v,w   lz      ly      lx
		self.vel  = np.zeros((3     , lr1[2], lr1[1], lr1[0]))
      #                    p       lz      ly      lx     
		self.pres = np.zeros((var[2], lr1[2], lr1[1], lr1[0]))
      #                    T       lz      ly      lx     
		self.temp = np.zeros((var[3], lr1[2], lr1[1], lr1[0]))
      #                    s_i     lz      ly      lx     
		self.scal = np.zeros((var[4], lr1[2], lr1[1], lr1[0]))

#==============================================================================
class exadata:
	"""
	    data
	    A class containing data for reading/writing binary simulation files
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
		self.elem   = [elem(var, lr1) for i in range(nel)]


#==============================================================================
def readdns(fname):
	"""
	    readdns(fname)
	    A function for reading binary data from the SIMSON binary format

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
	wdsz = 8
	realtype = 'd'
	#
	# identify endian encoding (and number passive scalar)
	etagb = infile.read(4)
	etagL = struct.unpack('<i', etagb)[0]
	etagB = struct.unpack('>i', etagb)[0]
	# 1644 = 44 + (2*8) * 100 means a maximum of 100 passive scalars
	if ((etagL >= 44) & (etagL <= 1644)):
		# print('Reading little-endian file\n')
		emode = '<'
		nscal = (etagL-44)/(2*wdsz)
	elif ((etagB >= 44) & (etagB <= 1644)):
		# print('Reading big-endian file\n')
		emode = '>'
		nscal = (etagB-44)/(2*wdsz)
	else:
		print('ERROR: could not initerpret endianness')
		return -3 
	#
	# Reynolds Number
	Re = infile.read(wdsz)
	Re = struct.unpack(emode+realtype,Re)[0]
	#
	# Poiseuille/Couette flag [deprecated]
	PouCou = infile.read(4)
	PouCou = struct.unpack(emode+'i',PouCou)[0]
	#
	# physical box size (here x and z only)
	boxsz = [0 for i in range(3)]
	dum = infile.read(2*wdsz)
	dum = list(struct.unpack(emode+2*realtype,dum))
	boxsz[0] = dum[0]; boxsz[2] = dum[1]
	#
	# time
	time = infile.read(wdsz)
	time = struct.unpack(emode+realtype,time)[0]
	#
	# dummy variable
	dum = infile.read(wdsz)
	#
	# passive scalar(s) (not tested)
	pr = np.zeros(nscal)
	m  = np.zeros(nscal)
	for i in range(nscal):
		pr[i] = infile.read(wdsz)
		pr[i] = struct.unpack(emode+realtype,pr[i])[0]
		m[i] = infile.read(wdsz)
		m[i] = struct.unpack(emode+realtype,m[i])[0]
	#
	# end-of-line
	eol = infile.read(8)	
	#
	# box size
	lr1 = infile.read(3*4)
	lr1 = list(struct.unpack(emode+3*'i',lr1))
	#
	# nfzsym (z-symmetry flag)
	nfzsym = infile.read(4)
	nfzsym = list(struct.unpack(emode+'i',nfzsym))[0]
	#
	# end-of-line
	eol = infile.read(8)
	#
	# compute total number of points per element
	npel = lr1[0] * lr1[1] * lr1[2]
	#
	# get number of pysical dimensions
	ndim = 2 + (lr1[2]>1)
	#
	# flow type
	fltype = infile.read(4)
	fltype = struct.unpack(emode+'i',fltype)[0]
	#
	# delta star (and boxsz along y-direction)
	dstar = infile.read(wdsz)
	dstar = struct.unpack(emode+realtype,dstar)[0]
	boxsz[1] = 2/dstar
	#
	# end-of-line
	eol = infile.read(8)
	#
	# flow-type dependent quantities
	#
	if (fltype == -1):
		rlam = infile.read(wdsz)
		rlam = struct.unpack(emode+realtype,rlam)[0]
		eol = infile.read(8)
	if (fltype == -2):
		rlam = infile.read(wdsz)
		rlam = struct.unpack(emode+realtype,rlam)[0]
		spanv = infile.read(wdsz)
		spanv = struct.unpack(emode+realtype,spanv)[0]
		eol = infile.read(8)
	if ((fltype == 4) or (fltype == 5)):
		bstart = infile.read(wdsz)
		bstart = struct.unpack(emode+realtype,bstart)[0]
		blength = infile.read(wdsz)
		blength = struct.unpack(emode+realtype,blength)[0]
		eol = infile.read(8)
	if ((fltype >=4) and (fltype <= 9)):
		bstart = infile.read(wdsz)
		bstart = struct.unpack(emode+realtype,bstart)[0]
		blength = infile.read(wdsz)
		blength = struct.unpack(emode+realtype,blength)[0]
		rlam = infile.read(wdsz)
		rlam = struct.unpack(emode+realtype,rlam)[0]
		spanv = infile.read(wdsz)
		spanv = struct.unpack(emode+realtype,spanv)[0]
		eol = infile.read(8)
	if (abs(fltype) == 20):
		gr = infile.read(nscal*wdsz)
		gr = struct.unpack(emode+nscal*realtype, gr)
		eol = infile.read(8)
	#
	# get variables
	var = [0 for i in range(5)]
	var[0] = ndim  # position
	var[1] = ndim  # velocity
	var[2] = 0     # pressure is not saved (SIMSON) 
	var[3] = 0     # temperature is treated like a scalar (SIMSON)
	var[4] = nscal # scalars
	#
	#---------------------------------------------------------------------------
	# READ DATA
	#---------------------------------------------------------------------------
	#
	# number of points
	npel = lr1[0]*lr1[1]*lr1[2]
	#
	# number of points per plane
	nppl = lr1[0]*lr1[2]
	#
	# reading buffer in fourier space
	fou =    np.zeros((lr1[2],lr1[1],lr1[0]/2+1)) + \
	      1j*np.zeros((lr1[2],lr1[1],lr1[0]/2+1))
	#
	# initialize data structure
	data = exadata(ndim, 1, lr1, var)
	data.time   = time
	data.wdsz   = wdsz
	if (emode == '<'):
		data.endian = 'little'
	elif (emode == '>'):
		data.endian = 'big'
	#
	# generate geometry
	# - x-direction
	dx = boxsz[0]/lr1[0]
	for ix in range(lr1[0]):
		data.elem[0].pos[0,:,:,ix] = dx * ix
	# - y-direction
	dy = np.arccos(-1.0)/(lr1[1]-1)
	for iy in range(lr1[1]):
		data.elem[0].pos[1,:,iy,:] = boxsz[1] * (1-np.cos(dy*iy))/2
	# - z-direction
	dz = boxsz[2]/lr1[2]
	for iz in range(lr1[2]):
		data.elem[0].pos[2,iz,:,:] = dz * (iz-lr1[2]/2)
	#
	# read velocity and transform in physical space
	for idim in range(3):
		for iz in range(lr1[2]):
			if (iz <= lr1[2]/2):
				izf = iz
			else:
				izf = lr1[2]/2 * 3 - (iz+1)
			for iy in range(lr1[1]):
				fi = infile.read(lr1[0]*wdsz)
				fi = list(struct.unpack(emode+lr1[0]*realtype, fi))
				ip = 0
				for ix in range(lr1[0]/2):
					fou[izf,iy,ix] = (fi[ip] + 1j*fi[ip+1]) * nppl * (-1)**idim
					ip += 2
				# end-of-line
				eol = infile.read(8)
		#
		# back to physical space
		data.elem[0].vel[idim,:,:,:] = np.fft.irfft2(fou,(lr1[0],lr1[2]),(2,0))
	#
	# read scalars and transform in physical space
	for ivar in range(var[4]):
		for iz in range(lr1[2]):
			if (iz <= lr1[2]/2):
				izf = iz
			else:
				izf = lr1[2]/2 * 3 - (iz+1)
			for iy in range(lr1[1]):
				fi = infile.read(lr1[0]*wdsz)
				fi = list(struct.unpack(emode+lr1[0]*realtype, fi))
				ip = 0
				for ix in range(lr1[0]/2):
					fou[izf,iy,ix] = (fi[ip] + 1j*fi[ip+1]) * nppl
					ip += 2
				# end-of-line
				eol = infile.read(8)
		#
		# back to physical space
		data.elem[0].scal[ivar,:,:,:] = np.fft.irfft2(fou,(lr1[0],lr1[2]),(2,0))
	#
	#---------------------------------------------------------------------------
	# CLOSE FILE 
	#---------------------------------------------------------------------------
	#
	# close file
	infile.close()
	#
	# output
	return data
