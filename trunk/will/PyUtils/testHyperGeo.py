import nose.tools
import PyHyperGeo
from decimal import *
import math

def testSimple():
	"""
	Test simple example HyperGeoPDF(10,100,50,40)
	"""
	val = PyHyperGeo.HyperGeoPDF(10,100,50,40)
	diff = abs(math.log(val, 10) - math.log(Decimal("3.5219e-5"),10))
	
	for tol in xrange(10):
		nose.tools.assert_true(diff < Decimal(str(10**(-tol))), 
				'Did not perform within %(tol)s : %(diff)s' % \
				{'tol': 10**-tol, 'diff': diff})

def testBadInputs():
	"""
	Tests for errors with improper inputs
	"""
	
	nose.tools.assert_raises(ValueError, PyHyperGeo.nchoosek, 10,100)
	
	nose.tools.assert_raises(ValueError, PyHyperGeo.HyperGeoPDF, 10,100,50,1)
	nose.tools.assert_raises(ValueError, PyHyperGeo.HyperGeoPDF, 10,100,50,600)
	nose.tools.assert_raises(ValueError, PyHyperGeo.HyperGeoPDF, 10,100,200,40)
	
	nose.tools.assert_raises(ValueError, PyHyperGeo.HyperGeoCDF, 10,100,50,1)
	nose.tools.assert_raises(ValueError, PyHyperGeo.HyperGeoCDF, 10,100,50,600)
	nose.tools.assert_raises(ValueError, PyHyperGeo.HyperGeoCDF, 10,100,200,40)
	
				
def CheckPDF(in_val):
	"""
	Testing Hard Vals HyperGeoPDF(k, N, m, n)
	"""
	
	val = PyHyperGeo.HyperGeoPDF(in_val[1], in_val[2],
										in_val[3], in_val[4])
	diff = Decimal(str(abs(math.log(val,10) - math.log(in_val[0],10))))
	
	for tol in xrange(10):
		nose.tools.assert_true(diff < Decimal(str(10**(-tol))), 
				'Did not perform within %(tol)s : %(diff)s' % \
				{'tol': 10**-tol, 'diff': diff})

def CheckCDF(in_val):
	"""
	Testing Hard Vals HyperGeoCDF(k, N, m, n)
	"""
	
	val = PyHyperGeo.HyperGeoCDF(in_val[1], in_val[2],
										in_val[3], in_val[4])
	diff = Decimal(str(abs(math.log(val,10) - math.log(in_val[0],10))))
	
	
	
	for tol in xrange(10):
		nose.tools.assert_true(diff < Decimal(str(10**(-tol))), 
				'Did not perform within %(tol)s : %(diff)s' % \
				{'tol': 10**-tol, 'diff': diff})

def testCheckPDF():
	"""
	Test HyperGeoPDF with overflow errors
	"""
	
	
	#known values are taken from matlab
	known_vals = [ (Decimal("1.075713318949161e-052"), 40 , 10700, 3382, 600),
					(Decimal("6.911079699056203e-045"), 50 , 10700, 3382, 600),
					(Decimal("4.823711570002099e-038"), 60 , 10700, 3382, 600),
					(Decimal("5.084558007964436e-032"), 70 , 10700, 3382, 600),
					(Decimal("1.022478843134504e-026"), 80 , 10700, 3382, 600),
					(Decimal("4.668619004638841e-022"), 90 , 10700, 3382, 600),
					(Decimal("5.536251349308503e-018"), 100 , 10700, 3382, 600),
					(Decimal("7.451516302018718e-003"), 170 , 10700, 3382, 600),
					(Decimal("2.490115000906640e-002"), 180 , 10700, 3382, 600),
					(Decimal("3.600131521508918e-002"), 190 , 10700, 3382, 600),
					(Decimal("2.304836196678693e-002"), 200 , 10700, 3382, 600),
					]
	for this_check in known_vals:
		yield CheckPDF, this_check
		
def testCheckCDF():
	"""
	Test HyperGeoCDF with overflow errors
	"""
	#known values are taken from matlab
	known_vals = [ (Decimal("1.255920284879654e-052"), 40 , 10700, 3382, 600),
					(Decimal("8.462746855220796e-045"), 50 , 10700, 3382, 600),
					(Decimal("6.950599166857216e-032"), 70 , 10700, 3382, 600),
					(Decimal("7.310526174195645e-022"), 90 , 10700, 3382, 600),
					(Decimal("9.399289534258735e-018"), 100 , 10700, 3382, 600),
					]
					
	for this_check in known_vals:
		yield CheckCDF, this_check

