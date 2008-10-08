from __future__ import with_statement
import itertools as IT
import HIVGenoTypeChecker as GenoTyping
import subprocess
import tempfile
import threading
import time
import datetime
import string
import os
import re
import types
import Queue
import string
import bisect
import PyVirus
import collections
import logging
from Bio import SeqIO
from Bio import SeqFeature as SeqFeatClass
#a hack but the only way I can get everything to work right
from Bio.SeqFeature import SeqFeature
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein, generic_nucleotide
from Bio.Blast import NCBIXML, NCBIStandalone
from GenomeDiagram import *
from reportlab.lib import colors
from AnnotUtils import *
from collections import defaultdict

from scipy.interpolate import interp1d
import numpy
import pylab

def constant_factory(value):
	return PatSeq(None, None, None, None)


class PatRecord():
	def __init__(self, DELTA_T, TYPE, VALUE):
		self.delta_t = DELTA_T
		self.type = TYPE
		self.val = VALUE
		
		
	def __eq__(self, RHS):
		comp1 = RHS.delta_t == self.delta_t
		comp2 = RHS.type == self.type
		comp3 = RHS.val == self.val
		return comp1 & comp2 & comp3
	
	def __gt__(self, RHS):
		return self.delta_t > RHS.delta_t
	
	def __lt__(self, RHS):
		return self.delta_t < RHS.delta_t
	
	def __hash__(self):
		return hash(str(self.delta_t) + str(self.type) + str(self.val))
		
class RXRegimine():
	def __init__(self, RX_NAMES, TIME_STARTED, TIME_STOPPED):
		self.rx_names = RX_NAMES
		self.time_started = TIME_STARTED
		self.time_stopped = TIME_STOPPED
	
	def __str__(self):
		out_str = str(self.rx_names) + ':' + str(self.time_started)
		return out_str




class PatSeq(PyVirus.ViralSeq):
	def __init__(self, PAT_ID, TEST_SEQ, SEQ_NAME, STUDY):
		#if SEQ_NAME is None then assume that TEST_SEQ is a SeqRecord
		self.seq_name = SEQ_NAME
		self.annotation = None
		self.feature_annot = []
		self.feature_annot_type = {}
		self.this_genome = None
		self.clin_time_line = []
		self.rx_time_line = []
		self.pat_id = None
		self.study = STUDY
		self.offset = None
		self.final_pos = None
		
		self.rna_tc = None
		self.cd4_tc = None
		
	def GenomeToBioPython(self):
		"""
		GenomeToBioPython
			Returns a generic SeqRecord in the biopython format for
		the whole genome sequence.
		"""
		return SeqRecord(Seq(self.my_sequence, generic_nucleotide),
						 id = str(self.study) + ':' + str(self.pat_id))
	
	def TranslateAll(self, REFBASE, WANTED_REF = None):
		"""
		TranslateAll
				Uses a BLASTx query to determine the which ORFs
		correspond to HIV-1 proteins.  It then stores the data in a
		dictionary keyed by "genename".
		"""
		blast_record = REFBASE.BLASTx(self.GenomeToBioPython())
		gene_dict = {}
		reg = re.compile('.*?\|(.*?):(\w{3}).*')
		aa_start = 0
		for this_check in blast_record.alignments:
			label = reg.match(this_check.title).groups()
			this_gene = str(label[1])
			this_ref = str(label[0])
			if (WANTED_REF != None) & (this_ref != WANTED_REF):
				#make sure we translate "in context" of a reference
				continue
			if this_check.hsps[0].score < 900:
				continue
			if not(gene_dict.has_key(this_gene)):
				ref_gene = REFBASE[this_ref].annotation[this_gene]
				
				aa_seq = str(this_check.hsps[0].query)
				start_pos = this_check.hsps[0].query_start
				end_pos = start_pos + this_check.length*3
				nt_seq = self.my_sequence[start_pos:end_pos]
				
				whole_genome_start_pos = start_pos + ref_gene.start#ref_gene.start
				whole_genome_end_pos = whole_genome_start_pos + this_check.length*3
				
				
				gene = Gene(this_gene, this_gene.upper(),
							nt_seq, aa_seq, start_pos, 
							end_pos)
				gene.GeneAnnot(None, aa_start, aa_start + len(aa_seq))
				
				aa_start += len(aa_seq)
				gene_dict[this_gene] = gene
			#end indent back
				# view_str = '%(seq)s:%(gene)s\n%(this)s\n%(ref)s' % \
								# {'seq':self.seq_name, 'gene':this_gene,
								# 'this':aa_seq,'ref':ref_gene.aa_seq}
				
				# logging.debug(view_str)
				
				
		self.annotation = gene_dict
		self.SetOffset(REFBASE[WANTED_REF])
	
	def DidTakeDrug(self, DRUG_SET):
		"""
		A utility function which will return True if the patient took the
		drug regimine described in DRUG_SET and false otherwise.
		"""
		
		if len(self.rx_time_line) == 0:
			print 'No drugs logged'
		
		for this_event in self.rx_time_line:
			if this_event.time_started.days == 0:
				print this_event
				if DRUG_SET <= this_event.rx_names:
					return True
		return False
	
	
	
	
	def SetOffset(self, REF_SEQ):
		"""
		Uses a REF_SEQ to create an offset for all features.  Currently uses
		a protein as an anchor and then creates an offset based on that.  
		Currently only works with SINGLE proteins in each sequence.
		"""
		
		offset = None
		last_gene = None
		for this_gene in REF_SEQ.annotation:
			if this_gene in self.annotation:
				this_offset = REF_SEQ.annotation[this_gene].start - \
								self.annotation[this_gene].start
				if last_gene == None:
					last_gene = this_gene
					offset = this_offset
				else:
					warn_str = '%(seq)s:too many genes mapped:%(this)s:%(last)s' \
						% {'seq':self.seq_name, 'this':this_gene, 'last':last_gene}
					logging.warning(warn_str)
					
					if this_offset < offset:
						last_gene = this_gene
						offset = this_offset
					
					
					
				
		
		if last_gene == None:
			warn_str = '%(seq)s:had no genes:PAT:%(my_genes)s:REF:%(ref_genes)s' \
						% {'seq':self.seq_name, 'my_genes':str(self.annotation.keys()), 
								'ref_genes':str(REF_SEQ.annotation.keys())}
			logging.warning(warn_str)
			self.offset = 0
		else:
			self.annotation[last_gene].start += offset
			self.annotation[last_gene].end += offset
			
			self.offset = offset
			self.final_pos = self.annotation[last_gene].end
		
		logging.debug(self.seq_name + ':offset:' + str(self.offset))
	
	def LogAnnotation(self, TYPE, POS, SEQ, HOM, NAME):
		"""
		LogAnnotation
			Logs various types of annotations about the sequence for various
			features such as Genes, homology, and ELMs ... use none if the data
			does not apply
		"""
		
		if self.offset == None:
			raise AttributeError, 'Offset is not set properly'
		
		if (type(POS) != types.TupleType):
			try:
				type(POS.type) == types.StringType
			except AttributeError:
				raise TypeError, 'If POS is not a tuple then it must be an ANNOT'
			
			new_feat = POS
			if POS.type != 'ELM':
				new_feat.start += self.offset
				new_feat.end += self.offset
		elif TYPE == 'Annot':
			new_feat = Annot(NAME, self.offset + POS[0], 
								self.offset + POS[1], None)
		elif TYPE == 'HomIsland':
			new_feat = HomIsland(SEQ, self.offset + POS[0], 
									self.offset + POS[1], HOM)
		elif TYPE == 'MIRNA':
			new_feat = HumanMiRNA(NAME, self.offset + POS[0], 
									self.offset + POS[1], HOM)
		elif TYPE == 'ELM':
			new_feat = ELM(NAME, self.offset + POS[0], 
							self.offset + POS[1], None)
		elif TYPE == 'TF':
			new_feat = TFSite(NAME, self.offset + POS[0], 
								self.offset + POS[1], HOM)
		else:
			raise KeyError
		
		if new_feat.start > 8000:
			logging.warning('Bad feat found:%(name)s:%(bad)s' % {'name':self.seq_name, 'bad':str(new_feat)})
		#else:
			#logging.warning('Good feat found:%(name)s:%(bad)s' % {'name':self.seq_name, 'bad':str(new_feat)})
		bisect.insort(self.feature_annot, new_feat)
		self.feature_annot_type[new_feat.type] = True
		
		
	def AnnotateClinical(self, DELTA_T, TYPE, VAL):
		"""
		A helper function which will add an annotation to the clinical 
		time-line.
		"""
		this_annot = PatRecord(DELTA_T, TYPE, VAL)
		bisect.insort(self.clin_time_line, this_annot)
		
		self.rna_tc = None
		self.cd4_tc = None
	
	def GetSequence(self):
		"""
		Looks through the self.clin_time_line to find the LAST sequence and
		sets self.my_sequence
		"""
		check_fun = lambda x: x.type == 'SEQ'
		seq_reqs = filter(check_fun, self.clin_time_line)
		if len(seq_reqs) != 0:
			self.my_sequence = seq_reqs[-1].val
		else:
			self.my_sequence = None
	
	def MakeNumpyTimeCourse(self):
		"""
		Creates a Numpy representation of the CD4 and RNA time-course for 
		faster interpolation.
		"""
		
		get_cd4 = lambda x: x.type == 'CD4'
		get_rna = lambda x: x.type == 'RNA'
		get_times = lambda x: x.delta_t.days//7
		get_nums = lambda x: x.val
		
		cd4_annots = filter(get_cd4, self.clin_time_line)
		rna_annots = filter(get_rna, self.clin_time_line)
		
		self.rna_tc = numpy.zeros((2,len(rna_annots)))
		self.cd4_tc = numpy.zeros((2,len(cd4_annots)))
		
		self.rna_tc[0,:] = map(get_times,rna_annots)
		self.rna_tc[1,:] = map(get_nums,rna_annots)
		
		self.cd4_tc[0,:] = map(get_times,cd4_annots)
		self.cd4_tc[1,:] = map(get_nums,cd4_annots)
		
		
	def CheckPatWindow(self, START, STOP):
		"""
		Checks to make sure the patient can evaluate over the START STOP
		window.
		"""
		if (self.rna_tc == None) | (self.cd4_tc == None):
			self.MakeNumpyTimeCourse()
		try:
			if (self.rna_tc[0,0] > START) | (self.cd4_tc[0,0] > START):
				return False
			elif (self.rna_tc[0,-1] < STOP) | (self.cd4_tc[0,-1] < STOP):
				return False
			else:
				return True
		except IndexError:
			return False
		
	
	
	def GetClinVal(self, TIME):
		"""
		A utility function which returns the CD4 and Viral-load of a patient
		at a specific time.  This handles the required interpolation.
		
		@param: TIME		The time of the desired sample.
		
		@returns:			A list [CD4, Viral-Load] for the desired sample.
		"""
		
		if (self.rna_tc == None) | (self.cd4_tc == None):
			self.MakeNumpyTimeCourse()
		
		#get data out if it is a deltatime object
		try:
			this_time = TIME.days//7
		except AttributeError:
			this_time = TIME
		
		try:
			rna_interp = interp1d(self.rna_tc[0,:], self.rna_tc[1,:], 
							copy = False)
			cd4_interp = interp1d(self.cd4_tc[0,:], self.cd4_tc[1,:], 
							copy = False)
			return numpy.append(cd4_interp(this_time), rna_interp(this_time))
		except ValueError:
			return numpy.array([numpy.nan, numpy.nan])
		except AssertionError:
			return numpy.array([numpy.nan, numpy.nan])
		
		
		
	def DetermineResponder(self, METHOD):
		"""
		Uses the clinical time-course to determine the 
		
		"""
		if METHOD == 'SD':
			#use the standard datum method
			val = self.GetClinVal(8)[1] // self.GetClinVal(0)[1]
			return val < 0.01
		elif METHOD == 'WND':
			times = [0,2,4,8,12,24]
			check_vals = map(lambda x: self.GetClinVal(x)[1], times)
			check_array = numpy.array(check_vals)
			check_array = numpy.diff(check_array) <= 0
			return numpy.sum(check_array) > 3
	
	def CheckFeatures(self, FEAT_LIST, FUDGE_FACTOR, NUM_CHECKS = 3, 
									CHECK_CUTOFF = 1):
		"""
		Checks a list of provided features and returns a boolean array 
		indicating whether that feature is present on this sequence.
		"""
		
		all_checks = numpy.zeros((len(FEAT_LIST), NUM_CHECKS))
		
		min_func = lambda x: x[0]
		
		for this_check_num in xrange(NUM_CHECKS):
			#do multiple checks with different permutations
			#it may be slightly order dependent
			check_feat = copy.deepcopy(self.feature_annot)
			numpy.random.shuffle(check_feat)
			
			#loop through each feature in FEAT_LIST
			for this_feat in xrange(len(FEAT_LIST)):
				poss_items = []
				#loop through each feature in the shuffled copy of this list
				
				if FEAT_LIST[this_feat].start + FUDGE_FACTOR < self.offset:
					continue
				if FEAT_LIST[this_feat].start - FUDGE_FACTOR > self.final_pos:
					continue
				
				for this_check in xrange(len(check_feat)):
					val = check_feat[this_check].FuzzyEquals(FEAT_LIST[this_feat],
													FUDGE_POS = FUDGE_FACTOR)
					#val == -1 refers to "wrong type or not within "
					if val != -1:
						poss_items.append((val, this_check))
				
				if len(poss_items) > 0:
					best_ind = min(poss_items, key = min_func)[1]
					all_checks[this_feat, this_check_num] = 1
					check_feat.pop(best_ind)
					
				
		avg_find = numpy.sum(all_checks, axis=1)
		
		found = avg_find > CHECK_CUTOFF
		
		return found
		
	def WriteFeatures(self, FEAT_LIST, FILE_HANDLE, FUDGE_FACTOR = 50, 
						NUM_CHECKS = 3, CHECK_CUTOFF = 1):
		"""
		Writes the feature list to a file-handle
		"""
		
		found_feats = self.CheckFeatures(FEAT_LIST, FUDGE_FACTOR, 
											NUM_CHECKS = NUM_CHECKS,
											CHECK_CUTOFF = CHECK_CUTOFF)
		
		trans_dict = {True: '1', False: '0'}
		found_feats = map(lambda x: trans_dict[x], found_feats.tolist())
		
		for this_feat in xrange(len(FEAT_LIST)):
			if FEAT_LIST[this_feat].start + FUDGE_FACTOR < self.offset:
				FILE_HANDLE.write('NaN\t')
			elif FEAT_LIST[this_feat].start - FUDGE_FACTOR > self.final_pos:
				FILE_HANDLE.write('NaN\t')
			else:
				FILE_HANDLE.write(found_feats[this_feat] + '\t')
		FILE_HANDLE.write('\n')
	
class PatBase(collections.defaultdict):
	"""
	PatBase is a subclass of collections.defaultdict.  The database is keyed
	by the PatID AS AN INTEGER!!!.
	The __missing__ function has been overloaded to create a new PatSeq 
	instance when an unknown ID is provided
	"""
	
	
	def __missing__(self, KEY):
		"""
		Initializes a PatSeq with the provided PatID if it could not find it 
		in the database
		"""
		self[KEY] = PatSeq(KEY, None, None, None)
		return self[KEY]
		
	
	
	def ReadDirec(self, DIREC):
		"""
		Reads the data in the directory and appends it into the self.pat_data 
		dictionary.
		"""
		
		def ProcessNUM(F_HANDLE, TYPE, STUDY):
			junk_line = F_HANDLE.next()
			for this_line in F_HANDLE:
				parts = string.split(this_line, '\t')
				pat_id = int(parts[0])
				this_date = datetime.timedelta(0,0,0,0,0,0, float(parts[1]))
				this_val = float(parts[2])
				
				this_pat = self[pat_id]
				
				this_pat.study = STUDY
				this_pat.AnnotateClinical(this_date, TYPE, this_val)
			
		def ProcessRX(F_HANDLE, STUDY):
			junk_line = F_HANDLE.next()
			parts = string.split(junk_line, '\t')
			
			DRUG_START = 5
			START_IND = 4
			STOP_IND = 5
			
			drug_names = parts[DRUG_START:]
			
			for this_line in F_HANDLE:
				parts = string.split(this_line, '\t')
				this_id = int(parts[0])
				
				this_start = datetime.timedelta(0,0,0,0,0,0,float(parts[START_IND]))
				this_stop = datetime.timedelta(0,0,0,0,0,0,float(parts[STOP_IND]))
				
				drugs = set()
				for this_part in zip(parts[DRUG_START:], drug_names):
					if this_part[0] == '1':
						drugs.add(this_part[1])
				
				this_drug_reg = RXRegimine(drugs, this_start, this_stop)
				
				this_pat = self[this_id]
				this_pat.rx_time_line.append(this_drug_reg)
				
		
		def ProcessFasta(F_HANDLE, STUDY):
			
			id_re = re.compile('PtID (\d*)')
			time_re = re.compile('Week (\d*)')
			fil_fun = lambda x: x in string.ascii_letters
			
			this_pat_id = None
			this_time = None
			this_seq = None
			for this_line in F_HANDLE:
				if this_line[0] == '>':
					this_pat_id = int(id_re.findall(this_line)[0])
					time_num = float(time_re.findall(this_line)[0])
					this_time = datetime.timedelta(0,0,0,0,0,0,time_num)
				else:
					this_seq = filter(fil_fun, this_line)
					
					this_pat = self[this_pat_id]
					this_pat.study = STUDY
					this_pat.pat_id = this_pat_id
					this_pat.seq_name = STUDY + ':' + str(this_pat_id)
					this_pat.AnnotateClinical(this_time, 'SEQ', this_seq.upper())
					
			
		
		present_files = os.listdir(DIREC)
		
		#file_fun = lambda x: x[0] != '.'
		#present_files = filter(file_fun, present_files)
		present_files.pop(0)
		
		for this_file in present_files:
			parts = string.split(this_file, '_')
			study = parts[0]
			f_type = parts[-1]
			
			with open(DIREC + this_file) as handle:
				if f_type == 'CD4.txt':
					ProcessNUM(handle, 'CD4', study)
				elif f_type == 'RNA.txt':
					ProcessNUM(handle, 'RNA', study)
				elif f_type == 'RX.txt':
					ProcessRX(handle, study)
				elif f_type == 'fasta.txt':
					ProcessFasta(handle, study)
				else:
					KeyError
		
		need_pop = []
		for this_pat in self:
			if self[this_pat].my_sequence == None:
				need_pop.append(this_pat)
		if len(need_pop) != 0:
			for this_pat in need_pop:
				self.pop(this_pat)
		
	def PatTimeCourseFig(self, FILE_NAME):
		"""
		Generates a set of scatter plots which plot the CD4 v Viral-Load
		of the patients through a time-course
		"""
		
		times = [0,2,4,8,12,24]
		#times = [0,8]
		
		base_axis = 230
		
		pat_list = self.keys()
		
		color_list = map(lambda x: x.DetermineResponder('WND'), self.values())
		
		logging.debug('checking: ' + str(color_list))
		
		resp_count = 0
		non_resp = 0
		for this_val in self.values():
			if this_val.DetermineResponder('WND'):
				resp_count += 1
			else:
				non_resp += 1
		logging.debug('Found %(r)d responders and %(nr)d non-resp' % \
						{'r':resp_count, 'nr':non_resp})
		
		
		
		text_dict = {'size':'smaller'}
		
		pylab.subplots_adjust(wspace=0.1, hspace=0.1)
		
		
		for this_time in xrange(1, len(times)+1):
			this_axes = pylab.subplot(base_axis + this_time)
			this_win = numpy.zeros((len(self),2))
			for i in xrange(len(pat_list)):
				if self[pat_list[i]].CheckPatWindow(times[0], times[-1]):
					this_win[i,:] = self[pat_list[i]].GetClinVal(times[this_time - 1])
			pylab.scatter(this_win[:,1], this_win[:,0], c=color_list)#, {'axes':this_axes})
			pylab.axis([0, 10, 0, 1500])
			
			
			if (this_time == 1) | (this_time == 4):
				pylab.ylabel('CD4 Count', text_dict)
			if this_time in [4,5,6]:
				pylab.xlabel('log(Viral Load)',text_dict)
			if this_time in [2,3]:
				pylab.xticks('')
				pylab.yticks('')
			if this_time in [1,2,3]:
				pylab.xticks('')
			if this_time in [5,6]:
				pylab.yticks('')
			
			
			
			pylab.title('%(time)d weeks' % {'time': times[this_time - 1]},text_dict)
		
		pylab.savefig(FILE_NAME)
		
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	