from __future__ import with_statement
import nose.tools
from testELM import SeqAns, FileSetup
from Bio import SeqIO
import itertools
import os
import re
import PyVirus
import PatUtils
import time
import datetime
import bisect
import PyBLAST
import AnnotUtils
import logging
import logging.handlers
import threading
import numpy
import string


def setUpModule():
	fmt='%(levelname)s \t %(threadName)s \t %(funcName)s \t %(asctime)s \t %(message)s'
	log_file = os.environ['PYTHONSCRATCH'] + 'testLog_KEEP.log'
	logging.basicConfig(level=logging.DEBUG,
						filename=log_file, 
						filemode='w',
						 format = fmt)

def testPatRecord():
	"""
	Test the PatRecord Class
	"""
	
	test_time = (1)
	
	test_rec = PatUtils.PatRecord(test_time, 'SEQ', 'ACTATACTA')
	
	nose.tools.assert_true(test_rec != None, 'Did not load PatRecord')
	nose.tools.assert_true(test_rec.delta_t != None, 'Did not save delta_t')
	nose.tools.assert_true(test_rec.type == 'SEQ', 'Did not save type')
	nose.tools.assert_true(test_rec.val == 'ACTATACTA', 'Did not save val')
	
def testPatRecordEQS():
	"""
	Test the PatRecord <, >, ==, bisect.insort and hash
	"""
	
	time1 = datetime.timedelta(1)
	time2 = datetime.timedelta(2)
	time3 = datetime.timedelta(3)
	
	test_rec1 = PatUtils.PatRecord(time1, 'SEQ', 'A1')
	test_rec2 = PatUtils.PatRecord(time2, 'SEQ', 'A2')
	test_rec3 = PatUtils.PatRecord(time3, 'SEQ', 'A3')
	
	test_rec1_2 = PatUtils.PatRecord(time1, 'SEQ', 'A1')
	
	nose.tools.assert_true(test_rec1 == test_rec1_2, 
							'__eq__ did not function properly')
	nose.tools.assert_true(hash(test_rec1) == hash(test_rec1_2), 
							'__hash__ did not function properly')
	nose.tools.assert_true(test_rec1 < test_rec2, 
							'__lt__ did not function properly')
	nose.tools.assert_true(test_rec3 > test_rec2, 
							'__gt__ did not function properly')
	
	test_list = []
	bisect.insort(test_list, test_rec1)
	bisect.insort(test_list, test_rec3)
	bisect.insort(test_list, test_rec2)
	bisect.insort(test_list, test_rec1_2)
	
	nose.tools.assert_true(test_list[0] == test_rec1, 'Wrong sorting')
	nose.tools.assert_true(test_list[1] == test_rec1_2, 'Wrong sorting')
	nose.tools.assert_true(test_list[2] == test_rec2, 'Wrong sorting')
	nose.tools.assert_true(test_list[3] == test_rec3, 'Wrong sorting')
	
	
def testPatSeq():
	"""
	Test the PatSeq object
	"""
	
	test_pat = PatUtils.PatSeq(12736, 'AACTCTATGATAGTA', 'test_seq', 'ACGT302')
	
	nose.tools.assert_true(test_pat != None, 'Could not load PatSeq')
	
def testPatBase():
	"""
	Test the PatBase and file reading.
	"""
	
	PAT_BASE_DIREC = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\stanfordDBs\\'
	
	
	test_base = PatUtils.PatBase()
	
	nose.tools.assert_true(test_base != None, 'Could not load PatBase')
	
	test_base.ReadDirec(PAT_BASE_DIREC)
	
	#a few PatIDs manually gathered for spot checking
	test_ids = [23424, 23430, 23434, 
				6527, 6543, 29405]
	test_study = ['ACTG320', 'ACTG320', 'ACTG320', 
					'ACTG302', 'ACTG302', 'ACTG302']
	
	for this_check in zip(test_ids, test_study):
		nose.tools.assert_true(test_base.has_key(this_check[0]), 
								'Could not find KEY: ' + str(this_check[0]))
		test_pat = test_base[this_check[0]]
		
		nose.tools.assert_true(test_pat.study == this_check[1],
			'Value had wrong STUDY ' + this_check[1] + ' : ' + test_pat.study)
	
def testMakeNumpyTimeCourse():
	"""
	Test the Numpy TimeCourse Conversion
	"""
	
	PAT_BASE_DIREC = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\stanfordDBs\\'
	
	test_base = PatUtils.PatBase()
	
	nose.tools.assert_true(test_base != None, 'Could not load PatBase')
	
	test_base.ReadDirec(PAT_BASE_DIREC)
	
	#a few PatIDs manually gathered for spot checking
	test_ids = [23424, 23430, 23434]
	this_cd4_0 = [178, 78, 17]
	
	for this_check in zip(test_ids, this_cd4_0):
		nose.tools.assert_true(test_base.has_key(this_check[0]), 
								'Could not find KEY: ' + str(this_check[0]))
		test_base[this_check[0]].MakeNumpyTimeCourse()
		
		test_cd4 = test_base[this_check[0]].cd4_tc
		
		nose.tools.assert_true(test_cd4[1,0] == this_check[1],
			'CD4 vals were not Numpy-ed properly got: %(g)d, expected: %(e)d'\
			% {'g':test_cd4[1,0], 'e':this_check[1]})
	
def testGetClinVal():
	"""
	Test the Clinical Value Interpolation
	"""
	
	PAT_BASE_DIREC = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\stanfordDBs\\'
	
	test_base = PatUtils.PatBase()
	
	nose.tools.assert_true(test_base != None, 'Could not load PatBase')
	
	test_base.ReadDirec(PAT_BASE_DIREC)
	
	#a few PatIDs manually gathered for spot checking
	test_ids = [23430, 23434]
	this_cd4_0 = [78, 17]
	this_cd4_4 = [140, 21]
	
	for this_check in zip(test_ids, this_cd4_0, this_cd4_4):
		nose.tools.assert_true(test_base.has_key(this_check[0]), 
								'Could not find KEY: ' + str(this_check[0]))
		
		#check to make sure interpolation does not happen outside 
		#of measured window
		test_ret = test_base[this_check[0]].GetClinVal(-100)
		nose.tools.assert_true(numpy.isnan(test_ret[0]), 
					'Did not return Nan for large Neg Time')
		
		test_ret = test_base[this_check[0]].GetClinVal(1100)
		nose.tools.assert_true(numpy.isnan(test_ret[0]), 
					'Did not return Nan for large Pos Time')
		
		test = test_base[this_check[0]].GetClinVal(2)
		
		within_win = (test[0] > this_check[1]) & (test[0] < this_check[2])
		
		
		nose.tools.assert_true(within_win, 
					'Interp value is not within window %(L)d < %(T)d < %(M)d' \
					% {'L':this_check[1], 'M':this_check[2], 'T':test[0]})
		
		
		
def testClinFig():
	"""
	Test Generating the Clinical Scatter Plot figure
	"""
	
	PAT_BASE_DIREC = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\stanfordDBs\\'
	
	test_base = PatUtils.PatBase()
	
	nose.tools.assert_true(test_base != None, 'Could not load PatBase')
	
	test_base.ReadDirec(PAT_BASE_DIREC)
	
	fig_name = os.environ['PYTHONSCRATCH'] + 'clin_figure.png'
	
	test_base.PatTimeCourseFig(fig_name)

def testDidTakeDrug():
	"""
	Test the Drug Checking
	"""
			
			
	
	PAT_BASE_DIREC = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\stanfordDBs\\'
	
	test_base = PatUtils.PatBase()
	
	nose.tools.assert_true(test_base != None, 'Could not load PatBase')
	
	test_base.ReadDirec(PAT_BASE_DIREC)
	
	present_files = os.listdir(PAT_BASE_DIREC)
	for this_file in present_files:
		parts = string.split(this_file, '_')
		study = parts[0]
		f_type = parts[-1]
		
		
		if  f_type == 'RX.txt':
			
			with open(PAT_BASE_DIREC + this_file) as handle:
				junk_line = handle.next()
				parts = string.split(junk_line, '\t')
				DRUG_START = 5
				START_IND = 4
				STOP_IND = 5
				
				drug_names = parts[DRUG_START:]
				
				
				
				for this_line in handle:
					
					parts = string.split(this_line, '\t')
					
					if float(parts[START_IND]) != 0:
						continue
					
					
					this_id = int(parts[0])
					for this_part in zip(parts[DRUG_START:], drug_names):
						super_set = set()
						counter = 0
						if this_part[0] == '1':
							correct = True
							super_set.add(this_part[1])
						else:
							correct = False
							counter += 1
							if counter == 10:
								counter = 0
							else:
								continue
								
						test_set = set()
						test_set.add(this_part[1])
						
						yield CheckDidTakeDrug, test_base, this_id, test_set, correct
					while len(super_set) > 2:
						yield CheckDidTakeDrug, test_base, this_id, super_set, True
						super_set.pop()
					
					
def CheckDidTakeDrug(PAT_BASE, PAT_ID, DRUG, CORRECT):
	
	bad_str = 'The PatID:%(pat_id)d ' % {'pat_id':PAT_ID}
	if CORRECT:
		bad_str += 'did not Log %(drug)s ' % {'drug':DRUG}
	else:
		bad_str += 'Logged %(drug)s ' % {'drug':DRUG}
		
	nose.tools.assert_true(PAT_BASE[PAT_ID].DidTakeDrug(DRUG) == CORRECT,
							bad_str)

							
def testGuessGene():
	"""
	Test whether the correct gene-fragment is translated
	"""
	
	base_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	dest_dir = os.environ['PYTHONSCRATCH']
	ref_base = PyVirus.RefBase(base_dir, dest_dir)
	
	for ref_seq in ref_base.ref_seqs:
		
		for this_gene in ref_seq.annotation:
			this_seq = ref_seq.annotation[this_gene].nt_seq
			logging.debug('Should be gene:' + this_gene)
			
			test_seq = PatUtils.PatSeq(1, None,
									'test','test_study')
			test_seq.my_sequence = this_seq
			test_seq.TranslateAll(ref_base, WANTED_REF = ref_seq.seq_name)
			
			found_genes = test_seq.annotation.items()
			found_genes = map(lambda x:x[0], found_genes)
			
			yield CheckCorrectGene, test_seq, found_genes, this_gene
			yield CheckNonFragment, test_seq, this_gene
			try:
				test_seq.annotation.pop(this_gene)
			except:
				pass
			for this_gene_new in test_seq.annotation:
				yield CheckIsFragment, test_seq, this_gene_new
			
			#yield CheckGuessGene, this_seq, ref_base, ref_seq.seq_name, this_gene

def CheckCorrectGene(test_seq, found_genes, CORRECT_GENE):
	nose.tools.assert_true(test_seq.annotation.has_key(CORRECT_GENE),
		'TranslateAll found %(found)s : needed %(correct)s' % \
		{'found': str(found_genes), 'correct':CORRECT_GENE})
	logging.debug('Found %(cor)s gene correctly' % {'cor':CORRECT_GENE})

def CheckNonFragment(test_seq, CORRECT_GENE):
	nose.tools.assert_true(test_seq.annotation[CORRECT_GENE].is_fragment == False,
		'TranslateAll incorrectly marked %(correct)s as a fragment!' % \
		{'correct':CORRECT_GENE})
	logging.debug('Found %(cor)s As NonFragment' % {'cor':CORRECT_GENE})

def CheckIsFragment(test_seq, this_gene):
	nose.tools.assert_true(test_seq.annotation[this_gene].is_fragment,
		'TranslateAll incorrectly marked %(correct)s as a whole gene!' % \
		{'correct':this_gene})
	logging.debug('Found %(cor)s as Fragment' % {'cor':this_gene})



							
							
#@nose.tools.with_setup(log_setup, None)
def testOffset():
	"""
	Test aligning fragments correctly
	"""
	#logging.getLogger('').setLevel(logging.CRITICAL)
	base_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	dest_dir = os.environ['PYTHONSCRATCH']
	ref_base = PyVirus.RefBase(base_dir, dest_dir)
	
	ELM_DICT = AnnotUtils.ELMParser()
	
	
	#since we don't want to deal with a time-course manually set the sequence!
	
	for ref_seq in ref_base.ref_seqs[0:3]:
		#ref_seq.FindELMs(ELM_DICT)
		ref_seq.FindTFSites()
		temp = map(lambda x: str(x), ref_seq.feature_annot)
		logging.debug('ref_seq tfs:' + str(temp))
		for this_gene in ref_seq.annotation:
			this_seq = ref_seq.annotation[this_gene].nt_seq
			logging.debug('should be gene:' + this_gene)
			
			test_seq = PatUtils.PatSeq(1, None,
									'test','test_study')
			test_seq.my_sequence = this_seq
			test_seq.TranslateAll(ref_base, WANTED_REF = str(ref_seq.seq_name))
			yield CheckOffset, ref_base, ref_seq, test_seq
			
		
		
	
def CheckOffset(REF_BASE, INPUT_REF, test_seq):
	
	
	
	test_seq.FindTFSites()
	
	
	if test_seq.offset == 0:
		logging.debug('BAD offset:' + str(test_seq.my_sequence))
	
	out_array = INPUT_REF.CheckFeatures(test_seq.feature_annot, 100)
	
	
	if not(out_array.all()):
		
		
		
		missing = []
		correct_pair = []
		
		for i in xrange(len(test_seq.feature_annot)):
			if out_array[i]==0:
				missing.append(str(test_seq.feature_annot[i]))
				
				
		nose.tools.assert_true(False, 
		'Found %(found)d of %(tot)d features: %(missing)s : %(found)s' % \
		{'found': numpy.sum(out_array), 'tot':numpy.size(out_array),
		'missing': str(missing)})
	
	
		
		
		
		
		
		
		
		
		