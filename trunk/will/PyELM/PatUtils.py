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
		
		self.rna_tc = None
		self.cd4_tc = None
		
		
		
	def __getattr__(self, NAME):
		if NAME == 'tested_subtype':
			self.DetSubtype()
			return self.tested_subtype
		elif NAME == 'my_sequence':
			self.GetSequence()
			return self.my_sequence
		else:
			raise AttributeError
		
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
		
		text_dict = {'size':'smaller'}
		
		pylab.subplots_adjust(wspace=0.1, hspace=0.1)
		
		
		for this_time in xrange(1, len(times)+1):
			this_axes = pylab.subplot(base_axis + this_time)
			this_win = numpy.zeros((len(self),2))
			for i in xrange(len(pat_list)):
				if self[pat_list[i]].CheckPatWindow(times[0], times[-1]):
					this_win[i,:] = self[pat_list[i]].GetClinVal(times[this_time - 1])
			print this_win
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
		
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	