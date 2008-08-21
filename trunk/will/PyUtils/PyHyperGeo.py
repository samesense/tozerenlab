"""	
can't calculate 1-HyperGeoPDF_safe(158,10700,3382,404) ... should be 4.7304e-4 but returns 0 or -inf

"""
from __future__ import with_statement
from decimal import *
import pickle
import sys
import itertools as IT

def nchoosek(n, k):
	"""
	returns the number of possible ways of choosing k items from a group of a 
	total of n items.
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
	
	k:	Num of genes in list with TF
	N:	Total num of genes in database
	m:	Num of genes with TF in database
	n:	Total Num of genes in list
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
	for i in xrange(1, k + 1):
		prob += HyperGeoPDF(i, N, m, n)
	
	return prob
	
def HyperGeoPDF(k, N, m, n):
	"""
	There is a shipment of N objects in which m are defective. HyperGeoCDF 
	returns the probability that EXACTLY k objects are defective in a sample 
	of n distinctive objects drawn from the shipment.
	
	TESTED SAFELY AGAINST OVERFLOW ERRORS!!!
	
	k:	Num of genes in list with TF
	N:	Total num of genes in database
	m:	Num of genes with TF in database
	n:	Total Num of genes in list
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