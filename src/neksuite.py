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
import exadata as exdat


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
	vars = header[-1].decode('utf-8')
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
	data = exdat.exadata(ndim, nel, lr1, var)
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
def readrea(fname):
	"""
	    readrea
	    A function for reading .rea files for nek5000

	    input variable:
	    fname : file name
	"""
	#
	try:
		infile = open(fname, 'r')
	except IOError as e:
		print('I/O error ({0}): {1}'.format(e.errno, e.strerror))
		#return -1
	#
	#---------------------------------------------------------------------------
	# READ HEADER (2 lines) + ndim + number of parameters
	#---------------------------------------------------------------------------
	#
	infile.readline()
	infile.readline()
	ndim = int(infile.readline().split()[0])
	npar = int(infile.readline().split()[0])
	#
	nface = 2*ndim
	#
	#---------------------------------------------------------------------------
	# READ parameters
	#---------------------------------------------------------------------------
	#
	param = np.zeros((npar,1))
	for ipar in range(npar):
		param[ipar] = float(infile.readline().split()[0])
	#
	#---------------------------------------------------------------------------
	# skip passive scalars
	#---------------------------------------------------------------------------
	#
	npscal = int(infile.readline().split()[0])
	for ipscal in range(npscal):
		infile.readline()
	#
	#---------------------------------------------------------------------------
	# skip logical switches
	#---------------------------------------------------------------------------
	#
	nswitch = int(infile.readline().split()[0])
	for iswitch in range(nswitch):
		infile.readline()
	#
	#---------------------------------------------------------------------------
	# skip XFAC,YFAC,XZERO,YZERO
	#---------------------------------------------------------------------------
	#
	infile.readline()
	#
	#---------------------------------------------------------------------------
	# READ MESH
	#---------------------------------------------------------------------------
	#
	infile.readline()
	nel = int(infile.readline().split()[0])
	#
	# initialize data structure
	lr1 = [2, 2, ndim-1]
	var = [ndim, 0, 0, 0, 0]
	#
	data = exdat.exadata(ndim, nel, lr1, var)
	#
	# read geometry
	data.lims.pos[:,0] =  float('inf')
	data.lims.pos[:,1] = -float('inf')
	for iel in range(nel):
		# skip element number and group
		infile.readline()
		for idim in range(var[0]-1): # if ndim == 3 do this twice
			for jdim in range(var[0]):
				fi = infile.readline().split()
				data.elem[iel].pos[jdim,idim,0,0] = float(fi[0])
				data.elem[iel].pos[jdim,idim,0,1] = float(fi[1])
				data.elem[iel].pos[jdim,idim,1,0] = float(fi[2])
				data.elem[iel].pos[jdim,idim,1,1] = float(fi[3])
	#
	#---------------------------------------------------------------------------
	# CURVED SIDE DATA
	#---------------------------------------------------------------------------
	#
	infile.readline()
	ncurved = int(infile.readline().split()[0])
	for icurved in range(ncurved):
		line = infile.readline()
		if (nel < 1e3):
			iedge = int(line[0:3])-1
			iel   = int(line[3:6])-1
			data.elem[iel].curv[iedge] = float(line[6:16])
		elif (nel < 1e6):
			iedge = int(line[0:2])-1
			iel   = int(line[2:8])-1
			data.elem[iel].curv[iedge] = float(line[8:18])
		else:
			iedge = int(line[0:2])-1
			iel   = int(line[2:12])-1
			data.elem[iel].curv[iedge] = float(line[12:22])
	#
	#---------------------------------------------------------------------------
	# BOUNDARY CONDITIONS (skip everything)
	#---------------------------------------------------------------------------
	#
	infile.readline()
	infile.readline()
	for iel in range(nel):
		for iface in range(nface):
			infile.readline()
	#
	#---------------------------------------------------------------------------
	# FORGET ABOUT WHAT FOLLOWS
	#---------------------------------------------------------------------------	
	#
	#
	# close file
	infile.close()
	#
	# output
	return data
