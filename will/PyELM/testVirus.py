from __future__ import with_statement
import nose.tools
from testELM import SeqAns, FileSetup
from Bio import SeqIO
import itertools
import os
import re
import PyVirus
import tempfile
import time
import PyBLAST
import AnnotUtils
import logging
import logging.handlers
import threading
import numpy
import copy

logger = logging.getLogger('')


def testLoading():
	"""
	Test Loading the Virus classes
	"""
	logger.info('Starting testLoading')
	possible_classes = [PyVirus.ViralSeq, PyVirus.BkgSeq]
	file_loc = os.environ['MYDOCPATH'] + 'PyELM\\'
	with open(file_loc + '50_seqs.fasta', mode = 'r') as file_handle:
		this_iter = itertools.izip(SeqIO.parse(file_handle, 'fasta'),
									iter(possible_classes),
									itertools.repeat(None, 5))
		for this_test in this_iter:
			yield CheckStringLoading, this_test[1], this_test[0]

		this_iter = itertools.izip(SeqIO.parse(file_handle, 'fasta'),
									iter(possible_classes),
									itertools.repeat(None, 5))
		for this_test in this_iter:
			yield CheckSeqRecordLoading, this_test[1], this_test[0]
        
def CheckStringLoading(INPUT_CLASS, INPUT_REC):
    genome = INPUT_CLASS(INPUT_REC.seq.tostring(), 'testFile')
    nose.tools.assert_not_equal(genome, None, 'Could not load from string')

def CheckSeqRecordLoading(INPUT_CLASS, INPUT_REC):
    genome = INPUT_CLASS(INPUT_REC, None)
    nose.tools.assert_not_equal(genome, None, 'Could not load from SeqRecord')

        

def testGenotyping_SLOW():
	"""
	Test Genotyping the Virus classes
	"""
	logger.info('Starting testGenotyping_SLOW')
	possible_classes = [PyVirus.ViralSeq, PyVirus.BkgSeq]
	file_loc = os.environ['MYDOCPATH'] + 'PyELM\\'
	file_loc = os.environ['MYDOCPATH'] + 'PyELM\\'
	with open(file_loc + '50_seqs.fasta', mode = 'r') as file_handle:	
		this_iter = itertools.izip(SeqIO.parse(file_handle, 'fasta'),
									iter(possible_classes))
		for this_test in this_iter:
			yield CheckGenotyping, this_test[1], this_test[0]
        
def CheckGenotyping(INPUT_CLASS, INPUT_REC):
    genome = INPUT_CLASS(INPUT_REC.seq.tostring(), 'testFile')
    genome.DetSubtype()
    nose.tools.assert_not_equal(genome, None)


def testRefLoading():
	"""
	Test Genbank parsing of Reference Sequences
	"""
	logger.info('Starting testRefLoading')
	base_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	filter_fun = lambda x: x[-3:] == '.gb'
	this_iter = itertools.izip(itertools.ifilter(filter_fun,
												iter(os.listdir(base_dir))),
								itertools.repeat(None,5))
	for this_file in this_iter:
		yield CheckRefLoading, base_dir+this_file[0]


def CheckRefLoading(INPUT_FILE):
    genome = PyVirus.RefSeq(INPUT_FILE)
    nose.tools.assert_not_equal(genome.annotation, None)

def testRefBaseLoading():
	"""
	Test RefBase loading
	"""
	logger.info('Starting testRefBaseLoading')
	base_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	dest_dir = os.environ['PYTHONSCRATCH']
	ref_base = PyVirus.RefBase(base_dir, dest_dir)
	nose.tools.assert_not_equal(ref_base, None, 'Could not Load Data')
	nose.tools.assert_not_equal(ref_base.ref_seqs, None, 'Data not present')


def testGenomeToBioPython():
	"""
	Test the whole genome conversion to BioPython SeqRecord
	"""
	logger.info('Starting testGenomeToBioPython')
	possible_classes = [PyVirus.ViralSeq, PyVirus.BkgSeq]
	file_loc = os.environ['MYDOCPATH'] + 'PyELM\\'
	with open(file_loc + '50_seqs.fasta', mode = 'r') as file_handle:
		this_iter = itertools.izip(SeqIO.parse(file_handle, 'fasta'),
									iter(possible_classes))
		for this_test in this_iter:
			yield CheckGenomeToBioPython, this_test[1], this_test[0]

   

def CheckGenomeToBioPython(INPUT_CLASS, INPUT_RECORD):

    genome = INPUT_CLASS(INPUT_RECORD.seq.tostring(), 'testseq')
    output = genome.GenomeToBioPython()
    nose.tools.assert_equal(output.id, 'testseq',
                            'Did not perform conversion properly')



def testBLASTContext():
	"""
	Test TempFile creation.
	"""
	logger.info('Starting testBLASTContext')
	with PyVirus.BLASTContext(os.environ['PYTHONSCRATCH'], 3) as file_names:
		for this_file in file_names:
			with open(this_file, mode = 'w') as handle:
				handle.write('testData')

	dir_list = os.listdir(os.environ['PYTHONSCRATCH'])
	for this_file in file_names:
		nose.tools.assert_false(this_file in dir_list,
								'Did not properly destroy files')
        

def testRefBaseBuilding():
    """
    Test Building the BLAST database
    """
    logger.info('Starting testRefBaseBuilding')
    correct_suffix_list = ['hr', 'in', 'sd', 'si', 'sq']
    correct_base_list = ['ref_aa.fasta.p', 'ref_nt.fasta.n']
    base_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
    dest_dir = os.environ['PYTHONSCRATCH']
    ref_base = PyVirus.RefBase(base_dir, dest_dir, BUILD = True)

    file_list = os.listdir(os.environ['PYTHONSCRATCH'])
    print file_list
    for this_base in correct_base_list:
        for this_suffix in correct_suffix_list:
            this_file = this_base + this_suffix
            nose.tools.assert_true(this_file in file_list,
                                   'Did not create:' + this_base + this_suffix)
    

def testTranslateAll():
	"""
	Test the translation of whole genome sequences
	"""
	logger.info('Starting testTranslateAll')
	base_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	dest_dir = os.environ['PYTHONSCRATCH']
	ref_base = PyVirus.RefBase(base_dir, dest_dir)
	possible_classes = [PyVirus.ViralSeq, PyVirus.BkgSeq]
	file_loc = os.environ['MYDOCPATH'] + 'PyELM\\'
	with open(file_loc + '50_seqs.fasta', mode = 'r') as file_handle:
		this_iter = itertools.izip(SeqIO.parse(file_handle, 'fasta'),
									iter(possible_classes))
		for this_test in this_iter:
			yield CheckTranslateAll, ref_base, this_test[1], this_test[0]
			yield CheckTrandlateAllCont, ref_base, this_test[1], this_test[0]

def CheckTranslateAll(REF_BASE, INPUT_CLASS, INPUT_RECORD):
    
    genome = INPUT_CLASS(INPUT_RECORD.seq.tostring(), 'testseq')
    genome.TranslateAll(REF_BASE)

    nose.tools.assert_true(genome.annotation.has_key('env'))   

	
def CheckTrandlateAllCont(REF_BASE, INPUT_CLASS, INPUT_RECORD):

	genome = INPUT_CLASS(INPUT_RECORD.seq.tostring(), 'testseq')
	wanted_ref = REF_BASE.ref_seqs[0].seq_name
	genome.TranslateAll(REF_BASE, WANTED_REF = wanted_ref)

	nose.tools.assert_true(genome.annotation.has_key('env'))   
	
	
def testGetSeqFeatures():
	"""
	Test the SeqFeature Output
	"""
	logger.info('Starting testGetSeqFeatures')
	base_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	dest_dir = os.environ['PYTHONSCRATCH']
	ref_base = PyVirus.RefBase(base_dir, dest_dir)
	possible_classes = [PyVirus.ViralSeq, PyVirus.BkgSeq]
	file_loc = os.environ['MYDOCPATH'] + 'PyELM\\'
	with open(file_loc + '50_seqs.fasta', mode = 'r') as file_handle:
		this_iter = itertools.izip(SeqIO.parse(file_handle, 'fasta'),
									iter(possible_classes))
		for this_test in this_iter:
			yield CheckSeqFeatures, ref_base, this_test[1], this_test[0]

def CheckSeqFeatures(REF_BASE, INPUT_CLASS, INPUT_RECORD):
	genome = INPUT_CLASS(INPUT_RECORD.seq.tostring(), 'testseq')
	genome.TranslateAll(REF_BASE)
	
	feature_list = genome.GetSeqFeatures()
	found_env = False
	for this_feat in feature_list:
		if this_feat.id == 'env':
			found_env = True
	
	nose.tools.assert_true(found_env, 'Could not find ENV in the feature list')   

def testAnnotClass():
	"""
	Test Annotation Classes
	"""
	logger.info('Starting testAnnotClass')
	base_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	dest_dir = os.environ['PYTHONSCRATCH']
	ref_base = PyVirus.RefBase(base_dir, dest_dir)
	
	pos_classes = [AnnotUtils.Annot, AnnotUtils.HumanMiRNA, 
					AnnotUtils.HomIsland, AnnotUtils.ELM,
					AnnotUtils.TFSite]
	
	
	for this_test in pos_classes:
		yield CheckAnnotClass, this_test
	
def CheckAnnotClass(TEST_CLASS):
	this_class = TEST_CLASS('test_name', 5, 10, 0.9)
	
	nose.tools.assert_true(this_class != None, 'Did not create class')
	
	this_feat = this_class.GetSeqFeature()
	print this_feat
	nose.tools.assert_true(this_feat != None, 'Did not convert properly')
	nose.tools.assert_true(this_feat.id == 'test_name', 
							'Did not convert values properly')
							

def testLogAnnotation():
	"""
	Test LogAnnotation
	"""
	logger.info('Starting testLogAnnotation')
	base_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	dest_dir = os.environ['PYTHONSCRATCH']
	ref_base = PyVirus.RefBase(base_dir, dest_dir)
	
	pos_classes = ['Annot', 'MIRNA', 'HomIsland']
	
	test_cases = itertools.izip(iter(ref_base.ref_seqs), 
								iter(pos_classes))
	for this_case in test_cases:
		yield CheckLogClass, this_case[0], this_case[1]

def CheckLogClass(THIS_SEQ, THIS_CLASS):
	THIS_SEQ.LogAnnotation(THIS_CLASS, (5, 10), 'aaaaa', 0.9, 'test_name')
	
	nose.tools.assert_true(len(THIS_SEQ.feature_annot) > 0, 
							'Did not annotate properly')
							
	this_feat = THIS_SEQ.feature_annot[0].GetSeqFeature()
	nose.tools.assert_true(this_feat.id != None, 
							'Did not convert values properly')

def testCheckFeatures():
	"""
	Test the Feature Checking
	"""
	logger.info('Starting testLogAnnotation')
	base_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	dest_dir = os.environ['PYTHONSCRATCH']
	
	ELM_DICT = AnnotUtils.ELMParser()
	
	ref_base = PyVirus.RefBase(base_dir, dest_dir)
	
	
	
	test_seq = ref_base.ref_seqs[0]
	test_seq.FindELMs(ELM_DICT)
	
	
	
	numpy.random.seed(22)
	
	#when intersected with itself it should return an array of ones
	
	for test_seq in ref_base.ref_seqs[1:]:
		test_seq.FindELMs(ELM_DICT)
		s = (1, len(test_seq.feature_annot))
		if len(test_seq.feature_annot) == 0:
			continue
		for i in range(5):
			bool_array = numpy.random.uniform(high = 10, size = s) > i
			yield CheckCheckFeatures, copy.deepcopy(test_seq), bool_array
	

def CheckCheckFeatures(TEST_SEQ, BOOL_ARRAY):
	
	testing_features = copy.deepcopy(TEST_SEQ.feature_annot)
	
	feats = numpy.where(BOOL_ARRAY==0)[1].tolist()
	feats.reverse()
	
	
	
	map(lambda x: TEST_SEQ.feature_annot.pop(x), feats)
	
	#TEST_SEQ.feature_annot = TEST_SEQ.feature_annot.pop(feats)
	
	output_array = TEST_SEQ.CheckFeatures(testing_features, 50)
	
	
	TP = numpy.sum(output_array[BOOL_ARRAY[0]])
	FP = numpy.sum(output_array[BOOL_ARRAY[0] == 0])
	TN = numpy.sum(output_array[BOOL_ARRAY[0] == 0] == 0)
	FN = numpy.sum(output_array[BOOL_ARRAY[0]] == 0)
	
	try:
		rec = float(TP)/float(TP + FN)
		pre = float(TP)/float(TP + FP)
	except ZeroDivisionError:
		raise ZeroDivisionError, 'TP:%(TP)d, FP:%(FP)d, TN:%(TN)d, FN:%(FN)d'%\
			{'TP': TP, 'FP': FP, 'TN': TN, 'FN': FN}
	
	f = 2 * float(pre * rec) / float(pre + rec)
	
	
	nose.tools.assert_true(f > 0.8, 
		'Found: F-Measure: %(F)f TP:%(TP)d, FP:%(FP)d, TN:%(TN)d, FN:%(FN)d'%\
		{'TP': TP, 'FP': FP, 'TN': TN, 'FN': FN, 'F':f})
	

def testMiRNA_SLOW():
	"""
	Test finding human miRNAs
	"""
	logger.info('Starting testMiRNA_SLOW')
	base_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	dest_dir = os.environ['PYTHONSCRATCH']
	ref_base = PyVirus.RefBase(base_dir, dest_dir)
	
	with open(dest_dir + 'RNAiCalibrations_KEEP.pkl') as handle:
			CALIB_DICT = pickle.load(handle)
	
	for this_calib in CALIB_DICT.keys()[10:]:
		junk = CALIB_DICT.pop(this_calib)
	
	this_seq = ref_base.ref_seqs[0]
	
	this_seq.HumanMiRNAsite(CALIB_DICT)
	
	num_feats = len(this_seq.feature_annot)
	
	nose.tools.assert_true(this_seq.feature_annot_type['MIRNA'], 
							'Did not set dictionary properly')
	nose.tools.assert_true(num_feats > 0, 'No Features annotated')
	
def testELM():
	"""
	Test finding ELMs
	"""
	logger.info('Starting testELM')
	base_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	dest_dir = os.environ['PYTHONSCRATCH']
	ref_base = PyVirus.RefBase(base_dir, dest_dir)
	
	ELM_DICT = AnnotUtils.ELMParser()
	
	test_iter = itertools.izip(iter(ref_base.ref_seqs), 
								itertools.repeat(ELM_DICT, 4))
	
	for this_test in test_iter:
		yield CheckELM, this_test[0], this_test[1]
		
def CheckELM(THIS_SEQ, ELM_DICT):
	THIS_SEQ.FindELMs(ELM_DICT)
	check_fun = lambda x: x.type == 'ELM'
	
	num_elms = filter(check_fun, THIS_SEQ.feature_annot)
	
	nose.tools.assert_true(num_elms > 0, 'Could not find any ELMs!')
	
def testTFFind_SLOW():
	"""
	Test finding TF binding sites
	"""
	logger.info('Starting testTFFind')
	base_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	dest_dir = os.environ['PYTHONSCRATCH']
	ref_base = PyVirus.RefBase(base_dir, dest_dir)
	
	test_iter = itertools.izip(iter(ref_base.ref_seqs), 
								itertools.repeat(None,4))
	
	for this_test in test_iter:
		yield CheckTF, this_test[0]
	
def CheckTF(THIS_SEQ):
	
	THIS_SEQ.FindTFSites()
	check_fun = lambda x: x.type == 'TFSite'
	
	num_tfs = filter(check_fun, THIS_SEQ.feature_annot)
	
	nose.tools.assert_true(num_tfs > 0, 'Could not find any TFs!')
	
def testGetSingleResult():
	"""
	Test the GetSingleResult of BLASTController
	"""
	logger.info('Starting testGetSingleResult')
	file_loc = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'
	with open(file_loc) as seq_handle:
		seq_iter = SeqIO.parse(seq_handle, 'fasta').next()
	
	b_control = PyBLAST.BLASTController(seq_iter, 'BLASTn', NUM_SEQS = 1)
	this_res = b_control.GetSingleResult()
	nose.tools.assert_true(this_res != None)

def testBLASTController_SLOW():
	"""
	Test the MultiThreaded BLAST controller
	"""
	logger.info('Starting testBLASTController_SLOW')
	file_loc = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'
	with open(file_loc) as seq_handle:
		seq_iter = SeqIO.parse(seq_handle, 'fasta')
		
		b_control = PyBLAST.BLASTController(seq_iter, 'BLASTn')
		before = time.time()
		b_control.start()
		t_diff = time.time() - before
		nose.tools.assert_true(t_diff < 20, 'Did not start in a new thread')
		
		counter = 0
		time.sleep(5)
		for this_rec in b_control.ResultGen():
			counter += 1
			nose.tools.assert_true(this_rec != None,
									'Returned bad value')
		nose.tools.assert_true(counter == 50, 
			'Did not return the proper number of sequences')




    
