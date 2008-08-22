from __future__ import with_statement
import os




class SeqDB():
	def __init__(self, CACHE_DIR, NAME = 'default.pkl'):
		self.cache_dir = CACHE_DIR
		self.name = NAME
	
	def IndexHIVSeqs(self):
		"""
		Reads the sequence files in the CACHE_DIR and indexes them for 
		searching.  The index is a defaultddict in which the KEY is a 
		n-mer sequence and the value is a list of tuples where:
			tup[0]		The sequence index
			tup[1]		The position of the n-mer
		"""
		pass
		
	def SearchHIVSeq(self, SEQ, NUM_MIS_MATCH):
		"""
		Searches the SEQ(s) within the database and their mismatches.
		
		@param SEQ:				A sequence (or list of sequences) to search
		@param NUM_MIS_MATCH	The number of mismatches to search
		@return:				A list which has the (n-mer found, seq_ind, pos)
		"""
		pass
		