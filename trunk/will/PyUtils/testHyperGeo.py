import nose.tools
import PyHyperGeo
from decimal import *

def testSimple():
	"""
	Test simple example HyperGeoPDF(10,100,50,40)
	"""
	val = PyHyperGeo.HyperGeoPDF(10,100,50,40)
	diff = abs(val - Decimal("3.5219e-5"))
	
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
	diff = abs(Decimal(str(val)) - in_val[0])
	
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
	diff = abs(Decimal(str(val)) - in_val[0])
	
	for tol in xrange(10):
		nose.tools.assert_true(diff < Decimal(str(10**(-tol))), 
				'Did not perform within %(tol)s : %(diff)s' % \
				{'tol': 10**-tol, 'diff': diff})

def testCheckPDF():
	"""
	Creates a generator for checking HyperGeoPDF with overflow errors
	"""
	
	
	
	known_vals = [ (Decimal("1.075713318949161e-052"), 40 , 10700, 3382, 600),
					(Decimal("6.911079699056203e-045"), 50 , 10700, 3382, 600),
					(Decimal("4.823711570002099e-038"), 60 , 10700, 3382, 600),
					(Decimal("5.084558007964436e-032"), 70 , 10700, 3382, 600),
					(Decimal("1.022478843134504e-026"), 80 , 10700, 3382, 600),
					(Decimal("4.668619004638841e-022"), 90 , 10700, 3382, 600),
					(Decimal("5.536251349308503e-018"), 100 , 10700, 3382, 600),
					(Decimal("1.896789735399537e-014"), 110 , 10700, 3382, 600),
					(Decimal("2.046831096279500e-011"), 120 , 10700, 3382, 600),
					(Decimal("7.469110706166772e-009"), 130 , 10700, 3382, 600),
					(Decimal("9.779853505596779e-007"), 140 , 10700, 3382, 600),
					(Decimal("4.830437366198409e-005"), 150 , 10700, 3382, 600),
					(Decimal("9.390576196102895e-004"), 160 , 10700, 3382, 600),
					(Decimal("7.451516302018718e-003"), 170 , 10700, 3382, 600),
					(Decimal("2.490115000906640e-002"), 180 , 10700, 3382, 600),
					(Decimal("3.600131521508918e-002"), 190 , 10700, 3382, 600),
					(Decimal("2.304836196678693e-002"), 200 , 10700, 3382, 600),
					]
	for this_check in known_vals:
		yield CheckPDF, this_check
		
def testCheckCDF():
	"""
	Creates a generator for checking HyperGeoCDF with overflow errors
	"""
		
	known_vals = [ (Decimal("1.255920284879654e-052"), 40 , 10700, 3382, 600),
					(Decimal("8.462746855220796e-045"), 50 , 10700, 3382, 600),
					(Decimal("6.224020458999711e-038"), 60 , 10700, 3382, 600),
					(Decimal("6.950599166857216e-032"), 70 , 10700, 3382, 600),
					(Decimal("1.490292674581295e-026"), 80 , 10700, 3382, 600),
					(Decimal("7.310526174195645e-022"), 90 , 10700, 3382, 600),
					(Decimal("9.399289534258735e-018"), 100 , 10700, 3382, 600),
					]
					
	for this_check in known_vals:
		yield CheckCDF, this_check

