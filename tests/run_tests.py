#------------------------------------------------------------------------------
# test nek scripts
#
def test_nek_scripts():
	import sys
	sys.path.append('./src/')
	import neksuite as ns
	#import numpy as np
	#import scipy.interpolate as spi
	#import matplotlib.pyplot as plt
	#from mayavi import mlab
	import time

	fname = './tests/nek/channel3D_0.f00001'

	# test file reading
	#
	ts = time.time()
	field = ns.readnek(fname)
	te = time.time()

	assert field.nel == 512

#------------------------------------------------------------------------------
# test simson scripts
#
def test_simson_scripts():
	import sys
	sys.path.append('./src/')
	import simsonsuite as ss
	#import numpy as np
	#import scipy.interpolate as spi
	#import matplotlib.pyplot as plt
	#from mayavi import mlab
	import time

	fname = './tests/simson/channel3D_t10000v.u'

	# test file reading
	#
	ts = time.time()
	field = ss.readdns(fname)
	te = time.time()

	assert field.nel == 1
	assert field.var == [3, 3, 0, 0, 0]
