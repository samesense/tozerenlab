"""
PyHyperGeo
	A set of functions for doing the hypergeometric test on large sets of 
	numbers without the worry of overflow/underflow/round-off errors.

"""
from decimal import *
import itertools as IT
import types


class memoized(object):
   """Decorator that caches a function's return value each time it is called.
   If called later with the same arguments, the cached value is returned, and
   not re-evaluated.
   """
   def __init__(self, func):
      self.func = func
      self.cache = {}
   def __call__(self, *args):
      try:
         return self.cache[args]
      except KeyError:
         self.cache[args] = value = self.func(*args)
         return value
      except TypeError:
         # uncachable -- for instance, passing a list as an argument.
         # Better to not cache than to blow up entirely.
         return self.func(*args)
   def __repr__(self):
      """Return the function's docstring."""
      return self.func.__doc__



def CheckList(SET_1, SET_2, BACKGROUND):
	"""
	Determines the p-value of list overlaps based on the HyperGeometric 
	distribution.
	
	@param:	SET_1			A set() of the elements which were 'picked' by an
							analysis method.
	
	@param:	SET_2			A set() of the elements which were 'picked' by a 
							seperate analysis method.
	
	@param:	BACKGROUND		A set() of the elements which COULD HAVE been
							picked by BOTH methods.
	
	@returns:				The p-value of the overlap based on the 
							hypergeometric distribution.
	
	Example:
		BACKGROUND:			All genes which have upstream promoters
		SET_1:				All genes which have the FOXD3 TF binding site
		SET_2:				All genes found by a SAM microarray analysis.
	
		returns:			The p-value of finding the overlap between the SAM
							analysis and the FOXD3 list.
	"""
	
	if type(SET_1) == types.ListType:
		raise TypeError, 'PICKED must be a set()'
	if type(SET_2) == types.ListType:
		raise TypeError, 'SET_2 must be a set()'
	if type(BACKGROUND) == types.ListType:
		raise TypeError, 'BACKGROUND must be a set()'
	
	if not(SET_1.issubset(BACKGROUND)):
		raise ValueError, 'SET_1 must be a subset of BACKGROUND'
	if not(SET_2.issubset(BACKGROUND)):
		raise ValueError, 'SET_2 must be a subset of BACKGROUND'
	
	N = len(BACKGROUND)
	k = len(SET_1 & SET_2)
	m = len(SET_2)
	n = len(SET_1)
	
	print (k - 1, N, m, n)
	
	return Decimal(1) - HyperGeoCDF(k - 1, N, m, n)
	
	
	


@memoized
def nchoosek(n, k):
	"""
	returns the number of possible ways of choosing k items from a group of a 
	total of n items.
	
	@param n:	The total number of items
	@param k:	The number of draws
	"""
	n = Decimal(n)
	k = Decimal(k)
	if 0 <= k <= n:  
		ntok = Decimal(1)
		ktok = Decimal(1)  
		for t in IT.imap(Decimal,xrange(1, min(k, n - k) + 1)):  
			ntok *= n  
			ktok *= t  
			n -= 1  
		return ntok / ktok
	else:  
		raise ValueError, 'n:%(n)s must be larger than k:%(k)s' %\
							{'n':n, 'k':k}
		
def HyperGeoCDF(k, N, m, n):
	"""
	There is a shipment of N objects in which m are defective. HyperGeoCDF 
	returns the probability that UPTO k objects are defective in a sample 
	of n distinctive objects drawn from the shipment.
	
	TESTED SAFELY AGAINST OVERFLOW/ROUND-OFF ERRORS!!!
	
	@param k:	Num of genes in list with TF
	@param N:	Total num of genes in database
	@param m:	Num of genes with TF in database
	@param n:	Total Num of genes in list
	"""
	
	if k > n:
		raise ValueError, 'n:%(n)s must be larger than k:%(k)s' %\
							{'n':n, 'k':k}
	elif n > N:
		raise ValueError, 'N:%(N)s must be larger than n:%(n)s' %\
							{'N':N, 'n':n}
	elif m > N:
		raise ValueError, 'N:%(N)s must be larger than m:%(m)s' %\
							{'N':N, 'm':m}
	
	
	prob = 0
	const_denom = nchoosek(N,n)
	n_ch_k_dict = {}
	
	for i in xrange(1, k + 1):
		if not(n_ch_k_dict.has_key((m,i))):
			n_ch_k_dict[(m,i)] = nchoosek(m,i)
		if not(n_ch_k_dict.has_key((N-m, n-i))):
			n_ch_k_dict[(N-m, n-i)] = nchoosek(N-m, n-i)
		
		prob += n_ch_k_dict[(m,i)]*n_ch_k_dict[(N-m, n-i)]/const_denom
	
	return prob

@memoized
def HyperGeoPDF(k, N, m, n):
	"""
	There is a shipment of N objects in which m are defective. HyperGeoCDF 
	returns the probability that EXACTLY k objects are defective in a sample 
	of n distinctive objects drawn from the shipment.
	
	TESTED SAFELY AGAINST OVERFLOW/ROUND-OFF ERRORS!!!
	
	@param k:	Num of genes in list with TF
	@param N:	Total num of genes in database
	@param m:	Num of genes with TF in database
	@param n:	Total Num of genes in list
	"""
	
	if k > n:
		raise ValueError, 'n:%(n)s must be larger than k:%(k)s' %\
							{'n':n, 'k':k}
	elif n > N:
		raise ValueError, 'N:%(N)s must be larger than n:%(n)s' %\
							{'N':N, 'n':n}
	elif m > N:
		raise ValueError, 'N:%(N)s must be larger than m:%(m)s' %\
							{'N':N, 'm':m}
	
	return nchoosek(m,k)*nchoosek(N-m,n-k)/nchoosek(N,n)
	