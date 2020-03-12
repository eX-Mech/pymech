#=============================================================================#
# exadata                                                                     #
#                                                                             #
#                                                                             #
#                                                                             #
# Authors: Jacopo Canton, Nicolo' Fabbiane, Guillaume Chauvat                 #
# Contacts: jacopo.canton(at)gmail.com                                        #
# Last edit: 2020-02-20                                                       #
#=============================================================================#
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
	def __init__(self, var, lr1, nbc):
		#                    x,y,z   lz      ly      lx
		self.pos  = np.zeros((3     , lr1[2], lr1[1], lr1[0]))
		#                    one per edge
		self.curv = np.zeros((12, 5))
		#             curvature type
		self.ccurv = ['' for i in range(12)]
		#                    u,v,w   lz      ly      lx
		self.vel  = np.zeros((3     , lr1[2], lr1[1], lr1[0]))
		#                    p       lz      ly      lx     
		self.pres = np.zeros((var[2], lr1[2], lr1[1], lr1[0]))
		#                    T       lz      ly      lx     
		self.temp = np.zeros((var[3], lr1[2], lr1[1], lr1[0]))
		#                    s_i     lz      ly      lx     
		self.scal = np.zeros((var[4], lr1[2], lr1[1], lr1[0]))
		#                    list of 8 parameters, one per face
		#                    one column for velocity, one for temperature, and one for each scalar
		self.bcs  = np.zeros((nbc, 6), dtype='U3, i4, i4, f8, f8, f8, f8, f8')


#==============================================================================
class exadata:
	"""
	    data
	    A class containing data for reading/writing binary simulation files
	"""
	def __init__(self, ndim, nel, lr1, var, nbc=0):
		self.ndim   = ndim
		self.nel    = nel
		self.ncurv  = []
		self.nbc    = nbc
		self.var    = var
		self.lr1    = lr1
		self.time   = []
		self.istep  = []
		self.wdsz   = []
		self.endian = []
		self.lims   = datalims(var)
		self.elem   = [elem(var, lr1, nbc) for i in range(nel)]

	def check_connectivity(self):
		dim = self.ndim
		err = False
		for (iel, el) in enumerate(self.elem):
			for ibc in range(self.nbc):
				for iface in range(2*dim):
					cbc = el.bcs[ibc, iface][0]
					if cbc == 'E' or cbc == 'P':
						connected_iel = int(el.bcs[ibc, iface][3])-1
						connected_face = int(el.bcs[ibc, iface][4])-1
						cbc1 = self.elem[connected_iel].bcs[ibc, connected_face][0]
						iel1 = int(self.elem[connected_iel].bcs[ibc, connected_face][3])-1
						iface1 = int(self.elem[connected_iel].bcs[ibc, connected_face][4])-1
						if iel1 < 0 or iel1 >= self.nel:
							err = True
							print("face {:} of element {:} is connected to face {:} of the nonexistent element {:}".format(iface, iel, connected_face, connected_iel))
						else:
							if cbc1 != cbc:
								err = True
								print("mismatched boundary conditions: face {:} of element {:} with condition {:} is connected to face {:} of element {:} with condition {:}".format(iface+1, iel+1, cbc, connected_face+1, connected_iel+1, cbc1))
							if iel1 != iel:
								err = True
								print("mismatched elements: face {:} of element {:} is connected to face {:} of element {:} but that face is connected to face {:} of element {:}".format(iface+1, iel+1, connected_face+1, connected_iel+1, iface1+1, iel1+1))
							if iface1 != iface:
								err = True
								print("mismatched faces: face {:} of element {:} is connected to face {:} of element {:} but that face is connected to face {:} of element {:}".format(iface+1, iel+1, connected_face+1, connected_iel+1, iface1+1, iel1+1))
		return not err

