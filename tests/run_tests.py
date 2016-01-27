#------------------------------------------------------------------------------
# test nek scripts
#
def test_nek_scripts():
	import sys
	sys.path.append('../src/')
	import neksuite as ns
	import numpy as np
	import scipy.interpolate as spi
	#import matplotlib.pyplot as plt
	#from mayavi import mlab
	import time

	fname = 'nek/channel3D_0.f00001'

	ts = time.time()
	field = ns.readnek(fname)
	te = time.time()

	if (field.nel == 18432):
		return 0
	else:
		return 1

#------------------------------------------------------------------------------
# test simson scripts
#
def test_simson_scripts():
	return 1
