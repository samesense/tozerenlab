from __future__ import with_statement
import nose.tools
import HIVDatabase
import PyVirus
import os
import time
import numpy
from Bio import SeqIO
import itertools
import re
import pickle
import AnnotUtils

def setupModule():
    dest_dir = os.environ['PYTHONSCRATCH']
    source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
    
    PyVirus.RefBase(source_dir, dest_dir, BUILD = True)



def testMappingRecord():
    """
    Test loading MappingRecord
    """
    dest_dir = os.environ['PYTHONSCRATCH']
    source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
    seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'
    seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

    with open(seq_file) as handle:
        seq_iter = itertools.imap(PyVirus.BkgSeq,
                                  SeqIO.parse(handle, 'fasta'),
                                  itertools.repeat(None, 5))
        for this_seq in seq_iter:
            yield CheckMappingRecord, this_seq


def CheckMappingRecord(INPUT_SEQ):
    this_rec = HIVDatabase.MappingRecord('testRef', INPUT_SEQ, 5)
    nose.tools.assert_equal(hash(this_rec), hash('testRef' + INPUT_SEQ.seq_name))

    
def testMappingBase():
    dest_dir = os.environ['PYTHONSCRATCH']
    source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'

    m_b = HIVDatabase.MappingBase(source_dir, dest_dir,
                                           'test_self', 5)

    nose.tools.assert_not_equal(m_b, None)
    

def testMappingBase():
    """
    Test loading MappingBase
    """
    dest_dir = os.environ['PYTHONSCRATCH']
    source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
    seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

    m_b = HIVDatabase.MappingBase(source_dir, 
                                           dest_dir, 'test_self')
    
    nose.tools.assert_not_equal(m_b, None)
    
def testAddtoShelf_SLOW():
	"""
	Test AddtoShelf
	"""
	dest_dir = os.environ['PYTHONSCRATCH']
	source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

	m_b = HIVDatabase.MappingBase(source_dir, dest_dir, 
											'test_self_KEEP')
	handle = open(seq_file)
	seq_iter = SeqIO.parse(handle, 'fasta')
	m_b.AddtoShelf(seq_iter)
	handle.close()
	
	for this_ref in m_b.ref_base:
		for this_test in m_b.test_names:
			this_key = this_ref.seq_name + this_test
			nose.tools.assert_true(m_b.my_map_shelf.has_key(this_key),
							'Shelf is missing an Item: ' + this_key)

def testSaveShelf():
	"""
	Test SaveShelf and autoloading
	"""

	dest_dir = os.environ['PYTHONSCRATCH']
	source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

	m_b = HIVDatabase.MappingBase(source_dir, dest_dir,
							   'test_self')
	POSSIBLE_KEY = 'AJ302647.1'
	with open(seq_file) as handle:
		seq_val = SeqIO.parse(handle, 'fasta').next()
		test_key = POSSIBLE_KEY + seq_val.id

	m_b.AddtoShelf([seq_val])

	del(m_b)

	m_b = HIVDatabase.MappingBase(source_dir, dest_dir,
							   'test_self')

	nose.tools.assert_true(m_b.my_map_shelf.has_key(test_key),
					'Shelf is missing an Item')

	for this_key in m_b.my_map_shelf:
		match_sum = numpy.sum(m_b.my_map_shelf[this_key].is_match)
		f_look_sum = numpy.sum(m_b.my_map_shelf[this_key].f_look)
		nose.tools.assert_true(match_sum > 0, 'There were no matchings found')
		nose.tools.assert_true(f_look_sum > 0, 'F_look was not properly calculated')

def testMappingSlicing():
	"""
	Test various slicing operations
	"""
	dest_dir = os.environ['PYTHONSCRATCH']
	source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

	m_b = HIVDatabase.MappingBase(source_dir, dest_dir,
									   'test_self')
	POSSIBLE_KEY = 'AJ302647.1'
	with open(seq_file) as handle:
		seq_val = SeqIO.parse(handle, 'fasta').next()
	test_key = POSSIBLE_KEY + seq_val.id

	m_b.AddtoShelf([seq_val])
	
	nose.tools.assert_true(m_b[test_key] != None,
							'Could not [] into MappingBase')
	test_mapping = m_b[test_key][0:5]
	nose.tools.assert_true(test_mapping != None,
							'Could not [:] into .is_match')

def testMappingBaseIter():
	"""
	Test the Iteration over MappingRecords
	"""
	
	dest_dir = os.environ['PYTHONSCRATCH']
	source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

	m_b = HIVDatabase.MappingBase(source_dir, dest_dir,
									   'test_self')
	POSSIBLE_KEY = 'AJ302647.1'
	with open(seq_file) as handle:
		seq_val = SeqIO.parse(handle, 'fasta').next()
	test_key = POSSIBLE_KEY + seq_val.id

	m_b.AddtoShelf([seq_val])

	for this_map in m_b.GetIter():
		nose.tools.assert_true(this_map != None, 
				'Could not generate individual mappings')
	
	for this_map in m_b.GetIter(WANTED_SUBTYPE = 'B'):
		nose.tools.assert_true(this_map != None, 
				'Could not generate individual mappings from a single subtype')
	
def testCheckSeqs():
	"""
	Test the Sequence Homology
	"""
	dest_dir = os.environ['PYTHONSCRATCH']
	source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

	m_b = HIVDatabase.MappingBase(source_dir, dest_dir,
									   'test_self')
	
	correct_seq = 'A'
	wrong_seq = 'QQQQQQ'
	
	nose.tools.assert_true(m_b.CheckSeqs(correct_seq) == 1,
							'Could not find "A" in all sequences.')
	nose.tools.assert_true(m_b.CheckSeqs(wrong_seq) == 0,
							'Found wrong_seq in a sequence')

def testMultiAnnot_SLOW():
	"""
	Test the multiThreaded Annotation
	"""
	
	dest_dir = os.environ['PYTHONSCRATCH']
	source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'
	
	with open(dest_dir + 'RNAiCalibrations_KEEP.pkl') as handle:
		CALIB_DICT = pickle.load(handle)
	for this_calib in CALIB_DICT.keys()[5:]:
		junk = CALIB_DICT.pop(this_calib)
	
	m_b = HIVDatabase.MappingBase(source_dir, dest_dir,
									   'test_self')
	
	ELM_DICT = AnnotUtils.ELMParser()
	wanted_ref = m_b.ref_base.ref_seqs[0]
	
	
	m_b.AnnotateBase(CALIB_DICT, ELM_DICT, True, wanted_ref)
	
	for this_seq in m_b.test_names:
		tf_val = m_b.my_test_shelf[this_seq].feature_annot_type['TF']
		mirna_val = m_b.my_test_shelf[this_seq].feature_annot_type['MIRNA']
		elm_val = m_b.my_test_shelf[this_seq].feature_annot_type['ELM']
		nose.tools.assert_true(tf_val, 'Did not process all TFs')
		nose.tools.assert_true(mirna_val, 'Did not process all MIRNAs')
		nose.tools.assert_true(elm_val, 'Did not process all ELMs')


def tearDownModule():
	checker = re.compile('.*KEEP.*')
	time.sleep(1)
	file_list = os.listdir(os.environ['PYTHONSCRATCH'])
	for this_file in file_list:
		if checker.match(this_file) == None:
			os.remove(os.environ['PYTHONSCRATCH'] + this_file)

    