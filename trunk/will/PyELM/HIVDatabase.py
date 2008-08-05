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
import string
import subprocess
import re
import string
from AnnotUtils import *
from GenomeDiagram import *
from reportlab.lib import colors

ref_source = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
dest_dir = "C:\\local_blast\\PyELMData\\"
bkg_file = os.environ['MYDOCPATH'] + 'PyELM\\500_seqs.fasta'

def enumerate(iterable):
    return itertools.izip(itertools.count(), iterable)
	
class MappingRecord():
	def __init__(self, REF_NAME, TEST_VIRAL, REF_LEN):
		self.ref_name = REF_NAME
		self.ref_len = REF_LEN

		self.test_name = TEST_VIRAL.seq_name
		self.test_seq_len = len(TEST_VIRAL.my_sequence)
		self.mapping = None
		self.is_match = None
		

	def __getitem__(self, KEY):
		return self.is_match[KEY]
	
	def __hash__(self):
		return hash(self.ref_name + self.test_name)

	def CalculateMapping(self, BLAST_ALIGNMENT):
		map_fun = lambda x: x != '-'
		fil_fun = lambda x: x < self.ref_len
		self.mapping = numpy.zeros((1,self.ref_len), dtype = numpy.uint16)
		self.is_match = numpy.zeros((1,self.ref_len), dtype = numpy.bool)
		filter_fun = lambda x: x[0] == '|'
		unzip = lambda x: x[1]
		for this_align in BLAST_ALIGNMENT.hsps:
			if this_align.query_start > this_align.query_end:
				continue
			if this_align.sbjct_start > this_align.query_end:
				continue
			if len(this_align.query) < 20:
				continue
			# query_inds = range(this_align.query_start, this_align.query_end)
			# subjct_inds = range(this_align.sbjct_start, this_align.sbjct_end)
			#print this_align.query
			#subjct_inds = range(this_align.sbjct_start, this_align.sbjct_end)
			#t_query_inds = range(this_align.query_start, this_align.query_end+1)
			query_inds = []
			subjct_inds = []
			counter=this_align.query_start
			#print this_align.query
			for this in xrange(len(this_align.query)):
				if this_align.query[this] != '-':
					counter+=1
				query_inds.append(counter)
			counter = this_align.sbjct_start
			for this in xrange(len(this_align.sbjct)):
				if this_align.sbjct[this] != '-':
					counter+=1
				subjct_inds.append(counter)
			subjct_inds = filter(fil_fun, subjct_inds)
			query_inds = subjct_inds[0:len(subjct_inds)]
			
			self.mapping[0,subjct_inds] = numpy.array(query_inds)
			matched_inds = filter(filter_fun, zip(this_align.match, query_inds))
			matched_inds = map(unzip, matched_inds)
			self.is_match[0,matched_inds] = 1
		
		last_val = 0
		for i in xrange(self.ref_len):
			if self.mapping[0,i] == 0:
				continue
			elif self.mapping[0,i] > last_val:
				last_val = self.mapping[0,i]
			else:
				self.mapping[0,i] = 0
				
	def CalculateRuns(self):
		"""
		CalculateRuns
			Find the number of consecutive MATCHES are present in the sequence
		"""
		p = numpy.nonzero(self.is_match)[1].tolist()
		self.r_look = numpy.zeros((1,self.ref_len), dtype = numpy.uint16)
		self.f_look = numpy.zeros((1,self.ref_len), dtype = numpy.uint16)
		for k,g in itertools.groupby(enumerate(p), lambda (i,x): i-x):
			inds = map(operator.itemgetter(1), g)
			#handle.write('inds \t' + str(inds))
			self.f_look[0,inds] = range(len(inds),0,-1)
			self.r_look[0,inds] = range(0, len(inds))
	
	def GetMapping(self, START_POS, END_POS):
		"""
		A utility function which converts the mapping position of the TEST_SEQ
		into the REF_SEQ
		"""
		try:
			#print self.mapping
			new_end = int(numpy.nonzero(self.mapping >= END_POS)[1][0])
		except IndexError:
			new_end = len(self.mapping)
		try:
			new_start = int(numpy.nonzero(self.mapping >= START_POS)[1][0])
		except IndexError:
			new_start = 0
		
		return new_start, new_end


class MappingBase():
	def __init__(self, REF_SOURCE, DEST_DIR, SHELF_NAME):
		self.ref_source = REF_SOURCE
		self.dest_dir = DEST_DIR
		self.ref_base = PyVirus.RefBase(REF_SOURCE, DEST_DIR)
		self.test_names = []
		self.my_test_shelf_name = DEST_DIR + SHELF_NAME + '_seqs.slf'
		self.my_map_shelf_name = DEST_DIR + SHELF_NAME + '_map.slf'
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
		
		if WANTED_REF == None:
			for this_ref in self.ref_base:
				if (WANTED_SUBTYPE == None) | (this_ref.tested_subtype == WANTED_SUBTYPE):
					for this_test in self.test_names:
						yield self.my_map_shelf[this_ref.seq_name+this_test]
	
		else:
			for this_test in self.test_names:
				yield self.my_map_shelf[WANTED_REF+this_test]



	def WindowedHomology(self, WANTED_REF):
		"""
		WindowedHomology
			Annotates a reference sequence with the windowed homology across 
			all sequences in the MAPPING_BASE based in the BLAST alignments.
		"""
		
		this_ref = self.ref_base.GetRefSeq(WANTED_REF)
		hom_vector = numpy.zeros((1, len(this_ref.my_sequence)), 
									'float')
		count = 0
		
		for this_mapping in self.GetIter(WANTED_REF = WANTED_REF):
			count += 1
			hom_vector += this_mapping.is_match
			
		hom_vector = hom_vector / float(count)
		
		val_iter = itertools.izip(range(len(hom_vector)),
									hom_vector.tolist())
		
		this_ref.global_hom = []
		for i in xrange(2,len(this_ref.my_sequence),2):
			this_ref.global_hom.append((i, hom_vector[0,i]))

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
			
			this_mapping = MappingRecord(ref_name, this_bkgseq, 
							len(self.ref_base.GetRefSeq(ref_name).my_sequence))
			this_mapping.CalculateMapping(this_align)
			this_mapping.CalculateRuns()
			
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
			
			
			
	def FindWindows(self, WANTED_REF, MIN_SIZE, MIN_HOM):
		
		num_seqs = float(len(self.test_names))
		
		this_ref = self.ref_base.GetRefSeq(WANTED_REF)
		start_pos = 0
		
		#pull all of the f_look arrays because getting each from the file 
		#is a time-intensive task
		multi_flook = numpy.zeros((len(self.test_names), 
									len(this_ref.my_sequence)),
									dtype = numpy.uint16)
		
		mapping_iter = self.GetIter(WANTED_REF)
		for i in xrange(len(self.test_names)):
			multi_flook[i,:] = mapping_iter.next().f_look
		
		while (start_pos + MIN_SIZE) < len(this_ref.my_sequence):
			
			#print win_size_array
			#print win_size_array.dtype
			vals = numpy.bincount(multi_flook[:,start_pos])
			#print vals
			
			win_size_norm = vals / num_seqs
			
			#win_size_prob[i] is the percentage of seqs with AT_LEAST i homology
			win_size_prob = numpy.flipud(numpy.cumsum(numpy.flipud(win_size_norm)))
			if len(win_size_prob) <= MIN_SIZE:
				start_pos += 1
				continue
			if win_size_prob[MIN_SIZE] < MIN_HOM:
				start_pos += 1
				continue
			else:
				
				for best_pos in xrange(MIN_SIZE, len(win_size_prob)):
					if win_size_prob[best_pos] < MIN_HOM:
						break
				if best_pos == 0:
					print 'found zero'
					raise KeyError
				this_seq = this_ref.my_sequence[start_pos:start_pos+best_pos]
				known_hom = self.CheckSeqs(this_seq)
				
				this_ref.LogAnnotation('HomIsland', (start_pos, 
										start_pos + best_pos),
										this_seq, known_hom, None)
				
				start_pos += best_pos
				
				
	def MakeMultiDiagram(self, WANTED_REF, FORCE = False):
		"""
		Uses GenomeDiagram to create a multi-genome figure.
		"""
		
		elm_filter = lambda x: x.type != 'ELM'
		
		this_ref = self.ref_base.GetRefSeq(WANTED_REF)
		
		this_gene_figure = GDDiagram(WANTED_REF)
		this_prot_figure = GDDiagram(WANTED_REF)
		
		dest_dir = os.environ['PYTHONSCRATCH']
		with open(dest_dir + 'RNAiCalibrations_KEEP.pkl') as handle:
			CALIB_DICT = pickle.load(handle)
		for this_calib in CALIB_DICT.keys()[10:]:
			junk = CALIB_DICT.pop(this_calib)
			
		ELM_DICT = ELMParser()
		
		gene_track = this_gene_figure.new_track(1, scale=0).new_set('feature')
		this_ref.WriteGenesToDiagram(gene_track)
		
		prot_tack = this_prot_figure.new_track(1, scale=0).new_set('feature')
		this_ref.WriteGenesToDiagram(prot_tack, PROT = True)
		
		#self.AnnotateBase(CALIB_DICT, ELM_DICT, FORCE, WANTED_REF)
		
		this_ref.HumanMiRNAsite(CALIB_DICT)
		this_ref.FindELMs(ELM_DICT)
		#this_ref.FindHomIslands(this_ref)
		this_ref.FindTFSites()
		
		#guess which feature
		#cached_seqs = self.my_test_shelf.values()
		all_features = []
		for this_seq in self.test_names:
			bkg_seq = self.my_test_shelf[this_seq]
			for this_annot in bkg_seq.feature_annot:
				if this_annot.CheckRange('ELM', 1000):
					all_features.append(this_annot)
				#if this_annot.CheckRange('HumanMiRNA', 7000):
					#all_features.append(this_annot)
		
		obj_vals = numpy.zeros((1,len(all_features)), dtype=numpy.uint16)
		
		current_min_val = 5000000
		current_annot = None
		pos_fun = lambda x: abs(x[0].start - x[1].start)
		for this_seq in self.test_names:
			bkg_seq = self.my_test_shelf[this_seq]
			
			for i in xrange(len(all_features)):
				this_annot = all_features[i]
				poss_annot = bkg_seq.FindEqAnnot(this_annot, 300)
				if poss_annot != None:
					obj_vals[0,i] += pos_fun((this_annot, poss_annot))
				else:
					obj_vals[0,i] += 1000
		
		best_inds = obj_vals.argsort()
		anchors = []
		#for i in xrange(1):
		anchors = copy.deepcopy(all_features[best_inds[0,0]])
		#anchors = copy.deepcopy(all_features[best_inds.tolist()])
		del(obj_vals)
		del(all_features)
		
		subtype_dict = {}
		
		for this_name in self.test_names:
			
			this_bkg_seq = self.my_test_shelf[this_name]
			this_subtype = this_bkg_seq.tested_subtype.split('|')[1]
			if subtype_dict.has_key(this_subtype):
				subtype_dict[this_subtype].append(this_name)
			else:
				subtype_dict[this_subtype] = [this_name]
			
			#save back into the shelve for faster calculation later
			self.my_test_shelf[this_name] = this_bkg_seq
		
		sorted_keys = subtype_dict.keys()
		sorted_keys.sort()
		sorted_keys = ['B', 'C']
		counter = 0
		for this_key in sorted_keys:
			for this_name in subtype_dict[this_key][0:min(len(subtype_dict[this_key]),300)]:
				this_bkg_seq = self.my_test_shelf[this_name]
				this_mapping = self.my_map_shelf[WANTED_REF + this_name]
				
				anchor_eq = this_bkg_seq.FindEqAnnot(anchors, 300, IS_LIST = False)
				
				#print 'here'
				if anchor_eq == None:
					continue
				
				this_gene_track = this_gene_figure.new_track(1, scale=0).new_set('feature')
				#this_prot_track = this_prot_figure.new_track(1, scale=0).new_set('feature')
				
				
				for this_annot in this_bkg_seq.feature_annot:
				
					# start_pos, end_pos = this_mapping.GetMapping(this_annot.start, 
																	# this_annot.end)
					# if (start_pos == 0) & (end_pos == 0):
						# continue
					# if start_pos == end_pos:
						# end_pos += 1
					
					new_annot = this_annot.MapMeToUS(anchor_eq, anchors)
					if new_annot != None:
						this_gene_track.add_feature(new_annot.GetSeqFeature(), 
											colour = new_annot.color)
					#print (str(this_annot), str(new_annot))
					#if this_annot.type == 'ELM':
					#	this_prot_track.add_feature(this_annot.GetSeqFeature_PROT(),
					#									colour = this_annot.color)
				
				counter += 1
				#raise KeyError
				#if counter > 500:
					#break
			#create a seperator
			#raise KeyError
			this_gene_figure.new_track(1, name = this_key).new_set('feature')
			this_prot_figure.new_track(1, name = this_key).new_set('feature')
		return (this_gene_figure, this_prot_figure)
		
		
	def AnnotateBase(self, MiRNA_DICT, ELM_DICT, FORCE, 
						WANTED_REF, WANTED_THREADS = 10):
		"""
		Makes multi-threaded calls to TFChecker.
		"""
		
		def AnnotWorker(SEQ, MI_DICT, E_DICT, FORCE_INPUT, REF_INPUT):
			
			SEQ.TranslateAll(self.ref_base, WANTED_REF = REF_INPUT)
			fil_fun = lambda x: x.type != 'ELM'
			SEQ.feature_annot = filter(fil_fun, SEQ.feature_annot)
			SEQ.FindELMs(E_DICT)
			
			
			if FORCE_INPUT or not(SEQ.feature_annot_type.get('MIRNA', False)):
				fil_fun = lambda x: x.type != 'MIRNA'
				SEQ
				SEQ.HumanMiRNAsite(MI_DICT)
			
			if FORCE_INPUT or not(SEQ.feature_annot_type.get('TFSite', False)):
				fil_fun = lambda x: x.type != 'TFSite'
				SEQ.FindTFSites()
			
			SEQ.DetSubtype()
			
			#make sure only one thread writes to the dictionary at a time
			#otherwise significant crap-age will happen
			with annot_lock:
				self.my_test_shelf[SEQ.seq_name] = SEQ
			worker_seph.release()
			
		annot_lock = threading.Lock()
		worker_seph = threading.Semaphore(WANTED_THREADS)
		all_threads = []
		
		for this_name in self.test_names:
			worker_seph.acquire()
			with annot_lock:
				#make sure only one thread writes to the dictionary at a time
				#otherwise significant crap-age will happen
				this_seq = self.my_test_shelf[this_name]
			this_thread = threading.Thread(target = AnnotWorker, 
								args = (this_seq, MiRNA_DICT, ELM_DICT, 
										FORCE, WANTED_REF))
			this_thread.start()
			all_threads.append(this_thread)
		#make sure all threads have closed before returning
		for this_thread in all_threads:
			this_thread.join()
		
		
		
		
		
		
		
		
		
		
		
		