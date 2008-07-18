from __future__ import with_statement
import PyVirus
import numpy
from Bio import SeqIO
import itertools
import copy
import re
import os
import threading
import copy
import logging
import cPickle as pickle
import types
import time
import shelve



ref_source = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
dest_dir = "C:\\local_blast\\PyELMData\\"
bkg_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

class MappingRecord():
	def __init__(self, REF_NAME, TEST_VIRAL):
		self.ref_name = REF_NAME

		self.test_name = TEST_VIRAL.seq_name
		self.test_seq_len = len(TEST_VIRAL.my_sequence)
		self.mapping = None
		self.is_match = None

	def __getitem__(self, KEY):
		return self.is_match[KEY]

	def __hash__(self):
		return hash(self.ref_name + self.test_name)

	def CalculateMapping(self, BLAST_ALIGNMENT):
		self.mapping = numpy.zeros((1,self.test_seq_len))
		self.is_match = numpy.zeros((1,self.test_seq_len), bool)
		filter_fun = lambda x: x[0] == '|'
		unzip = lambda x: x[1]
		for this_align in BLAST_ALIGNMENT.hsps:
			if this_align.query_start > this_align.query_end:
				continue
			if this_align.sbjct_start > this_align.query_end:
				continue	
			query_inds = range(this_align.query_start, this_align.query_end)
			subjct_inds = range(this_align.sbjct_start, this_align.sbjct_end)
			self.mapping[0,query_inds] = numpy.array(subjct_inds)

			matched_inds = filter(filter_fun, zip(this_align.match, query_inds))
			matched_inds = map(unzip, matched_inds)
			self.is_match[0,matched_inds] = 1





class MappingBase():
	def __init__(self, REF_SOURCE, DEST_DIR, SHELF_NAME):
		self.ref_source = REF_SOURCE
		self.dest_dir = DEST_DIR
		self.ref_base = PyVirus.RefBase(REF_SOURCE, DEST_DIR)
		self.test_names = []
		self.my_test_shelf_name = DEST_DIR + SHELF_NAME + '.slf'
		self.my_map_shelf_name = DEST_DIR + SHELF_NAME + '.seq'
		
		self.my_map_shelf = shelve.DbfilenameShelf(self.my_map_shelf_name, 
													protocol = 2)
		self.my_test_shelf = shelve.DbfilenameShelf(self.my_test_shelf_name, 
													protocol = 2)
		self.test_names = self.my_test_shelf.keys()

	
	def __getitem__(self, key):
		return self.my_map_shelf[key]

	def __setitem__(self, key, value):
		self.my_map_shelf[key] = value


	def BuildRefBase(self):
		"""
		BuildRefBase
		Build the reference BLAST databases.
		"""
		self.ref_base.BuildDatabase()

	def AddtoShelf(self, INPUT_SEQRECORD_ITER, MULTI_THREAD = False):
		"""
		AddtoShelf
			Blasts the background sequences against the RefBase and creates
		shelf instances of the resulting alignments.  Input must be a LIST,
		ITER or GENERATOR
		"""

		blast_seq_iter, map_seq_iter = itertools.tee(INPUT_SEQRECORD_ITER, 2)
		for this_blast in self.ref_base.BLASTn(blast_seq_iter, 
												NUM_THREADS = 4):
			
			for this_mapping in self.MapToRefs(this_blast, map_seq_iter.next()):
				volume_name = this_mapping.ref_name
				volume_name += this_mapping.test_name
				
				self.my_map_shelf[volume_name] = this_mapping

	def GetIter(self, WANTED_SUBTYPE = None):
		"""
		GetIter
			Returns an iterator over the mapping records in the database.  Can
			specifiy specific subtypes if desired.
		"""
		
		for this_ref in self.ref_base:
			if (WANTED_SUBTYPE == None) | (this_ref.tested_subtype == WANTED_SUBTYPE):
				for this_test in self.test_names:
					yield self.my_map_shelf[this_ref.seq_name+this_test]
	

	def MapToRefs(self, BLAST_RECORDS, SEQ_RECORD):
		"""
			MapToRefs
			Use a sequence of BLAST_RECORDs to create a sequence of
		MappingRecord.
		"""
		#reg_exp gets the Genbank ID from the header
		checker = re.compile('.*?\|(\w{1,2}\d*\.\d*).*')
		
		this_bkgseq = PyVirus.BkgSeq(SEQ_RECORD, None)
		if not(this_bkgseq.seq_name in self.test_names):
			self.test_names.append(this_bkgseq.seq_name)
			self.my_test_shelf[this_bkgseq.seq_name] = this_bkgseq
		
		for this_align in BLAST_RECORDS.alignments:
			temp_name = checker.match(this_align.title)
			ref_name = str(temp_name.groups()[0])
			
			this_mapping = MappingRecord(ref_name, this_bkgseq)
			this_mapping.CalculateMapping(this_align)
			
			yield this_mapping

