"""	
can't calculate 1-HyperGeoPDF_safe(158,10700,3382,404) ... should be 4.7304e-4 but returns 0 or -inf

"""

from __future__ import with_statement
import string
import os
import re
import threading
import logging
import PyMozilla
import pickle
import sys
import itertools as IT

def nchoosek(n, k):
	"""
	returns the number of possible ways of choosing k items from a group of a 
	total of n items.
	"""
	if 0 <= k <= n:  
		ntok = 1  
		ktok = 1  
		for t in xrange(1, min(k, n - k) + 1):  
			ntok *= n  
			ktok *= t  
			n -= 1  
		return ntok // ktok  
	else:  
		return 0
		
def nchoosek_gen(n,k):
	"""
	Returns a generator object which returns sucessive values of the 
	combinitoric data
	"""
	n_iter = IT.chain(xrange(n,1,-1), IT.repeat(1))
	k_iter = IT.chain(xrange(k,1,-1), IT.repeat(1))
	n_minus_k_iter = IT.chain(xrange(n-k,1,-1), IT.repeat(1))
	for vals in IT.izip(n_iter, k_iter, n_minus_k_iter):
		val = float(vals[0])/(float(vals[1])*float(vals[2]))
		yield val

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
	
	prob = 0
	for i in xrange(1, k + 1):
		prob += HyperGeoPDF_safe(i, N, m, n)
	
	return prob
	
def HyperGeoPDF_safe(k, N, m, n):
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

	m_ch_k_gen = nchoosek_gen(m,k)
	N_min_m_ch_n_min_k_gen = nchoosek_gen(N-m,n-k)
	N_ch_n_gen = nchoosek_gen(N,n)
	final_val = 1

	while True:
		if final_val > sys.maxint*0.1:
			final_val /= N_ch_n_gen.next()
		if final_val < 1.0e-400:
			num1 = m_ch_k_gen.next()
			num2 = N_min_m_ch_n_min_k_gen.next()
			final_val *= num1*num2
		else:
			num1 = m_ch_k_gen.next()
			num2 = N_min_m_ch_n_min_k_gen.next()
			denom = N_ch_n_gen.next()
			final_val *= num1*num2/denom
			if (num1 == 1) & (num2 == 1) & (denom == 1):
				break
	return final_val
