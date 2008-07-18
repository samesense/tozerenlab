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
import operator


ref_source = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
dest_dir = "C:\\local_blast\\PyELMData\\"
bkg_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

def enumerate(iterable):
    return itertools.izip(itertools.count(), iterable)

class MappingRecord():
	def __init__(self, REF_NAME, TEST_VIRAL):
		self.ref_name = REF_NAME

		self.test_name = TEST_VIRAL.seq_name
		self.test_seq_len = len(TEST_VIRAL.my_sequence)
		self.mapping = None
		self.is_match = None
		

	def __getitem__(self, KEY):
		return self.is_match[KEY]
	
	def __getattr__(self, NAME):
		if NAME == 'f_look':
			self.CalculateRuns()
			return self.f_look
		elif NAME == 'r_look':
			self.CalculateRuns()
			return self.r_look
		else:
			raise AttributeError
	
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
	
	def CalculateRuns(self):
		"""
		CalculateRuns
			Find the number of consecutive MATCHES are present in the sequence
		"""
		p = test_array.nonzero()[0].tolist()
		self.r_look = nump.zeros((1,self.test_seq_len))
		self.f_look = nump.zeros((1,self.test_seq_len))
		
		for k,g in itertools.groupby(enumerate(p), lambda (i,x): i-x):
			inds = map(operator.itemgetter(1), g)
			self.f_look[inds] = range(len(inds),0,-1)
			self.r_look[inds] = range(0, len(inds))
		
		




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

	def GetIter(self, WANTED_REF = None, WANTED_SUBTYPE = None):
		"""
		GetIter
			Returns an iterator over the mapping records in the database.  Can
			specifiy specific subtypes if desired.
		"""
		
		if WANTED_REF == None
			for this_ref in self.ref_base:
				if (WANTED_SUBTYPE == None) | (this_ref.tested_subtype == WANTED_SUBTYPE):
					for this_test in self.test_names:
						yield self.my_map_shelf[this_ref.seq_name+this_test]
	
		else:
			for this_test in self.test_names:
				yield self.my_map_shelf[WANTED_REF+this_test]

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

	def CheckSeqs(self, INPUT_SEQ):
		"""
		CheckSeqs
			Checks all test sequences in the database and returns a homology 
			percentage between 0.0 and 1.0
		"""
		counter = 0
		
		for this_seq in self.my_test_shelf:
			if self.my_test_shelf[this_seq].CheckSeq(INPUT_SEQ):
				counter += 1
		return float(counter) / float(len(self.test_names))
			
			
			
	def FindWindows(self, MIN_SIZE, MIN_HOM):
		
		
		initial_window = 1
		f_look_fun = opertator.attrgetter('f_look')
		pos_fun = lambda (x,y): x[y]
		for this_ref in self.ref_base:
			this_window_size = initial_window
			start_pos = 0
			poss_win = []
			while (start_pos + window_size) < len(this_ref.my_sequence):
				f_look_iter = itertools.imap(f_look_fun, 
								self.GetIter(WANTED_REF = this_ref.name))
				win_size_iter = itertools.imap(pos_fun, f_look_iter, start_pos)
				
				win_size_array = numpy.array(list(win_size_iter))
				
				win_size_count = numpy.bincount(win_size_array)
				
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	