"""
PyHyperGeo
	A set of functions for doing the hypergeometric test on large sets of 
	numbers without the worry of overflow/underflow/round-off errors.

"""
from decimal import *
import itertools as IT



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
	