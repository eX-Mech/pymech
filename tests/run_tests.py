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

	assert field.endian == 'little'
	assert field.istep  == 10
	assert field.lr1    == [8, 8, 8]
	assert field.ndim   == 3
	assert field.nel    == 512
	assert field.var    == [3, 3, 1, 0, 0]
	assert field.wdsz   == 4
	assert (field.time - 0.2) < 1e-3


	fname = './tests/nek/2D_section_R360.rea'

	# test .rea reading
	#
	ts = time.time()
	field = ns.readrea(fname)
	te = time.time()

	assert field.lr1  == [2, 2, 1]
	assert field.ndim == 2
	assert field.nel  == 1248
	assert (field.elem[0].pos[0][0][0][0] - 0.048383219999999998 ) < 1e-3


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

	assert field.endian == 'little'
	assert field.istep  == []
	assert field.lr1    == [48, 65, 48]
	assert field.ndim   == 3
	assert field.nel    == 1
	assert field.var    == [3, 3, 0, 0, 0]
	assert field.wdsz   == 8
	assert (field.time - 10000.439742009798) < 1e-3
