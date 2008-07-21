from __future__ import with_statement
import nose.tools
import HIVDatabase
import PyVirus
import os
import time
import numpy
from Bio import SeqIO
import itertools

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
    this_rec = HIVDatabase.MappingRecord('testRef', INPUT_SEQ)
    nose.tools.assert_equal(hash(this_rec), hash('testRef' + INPUT_SEQ.seq_name))

    
def testMappingBase():
    dest_dir = os.environ['PYTHONSCRATCH']
    source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'

    mapping_base = HIVDatabase.MappingBase(source_dir, dest_dir,
                                           'test_self')

    nose.tools.assert_not_equal(mapping_base, None)
    

def testMappingBase():
    """
    Test loading MappingBase
    """
    dest_dir = os.environ['PYTHONSCRATCH']
    source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
    seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

    mapping_base = HIVDatabase.MappingBase(source_dir, 
                                           dest_dir, 'test_self')
    
    nose.tools.assert_not_equal(mapping_base, None)
    
def testAddtoShelf_SLOW():
	"""
	Test AddtoShelf
	"""
	dest_dir = os.environ['PYTHONSCRATCH']
	source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

	mapping_base = HIVDatabase.MappingBase(source_dir, dest_dir, 'test_self')
	POSSIBLE_KEY = 'AJ302647.12_bkg_data_1'
	handle = open(seq_file)
	seq_iter = SeqIO.parse(handle, 'fasta')
	mapping_base.AddtoShelf(seq_iter)
	handle.close()

	nose.tools.assert_true(mapping_base.my_map_shelf.has_key(POSSIBLE_KEY), 'Shelf is missing an Item')

def testSaveShelf():
    """
    Test SaveShelf and autoloading
    """

    dest_dir = os.environ['PYTHONSCRATCH']
    source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
    seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

    mapping_base = HIVDatabase.MappingBase(source_dir, dest_dir,
                                           'test_self')
    POSSIBLE_KEY = 'AJ302647.1'
    with open(seq_file) as handle:
        seq_val = SeqIO.parse(handle, 'fasta').next()
	test_key = POSSIBLE_KEY + seq_val.id

	mapping_base.AddtoShelf([seq_val])

    del(mapping_base)

    mapping_base = HIVDatabase.MappingBase(source_dir, dest_dir,
                                           'test_self')
    
    nose.tools.assert_true(mapping_base.my_map_shelf.has_key(test_key),
                                'Shelf is missing an Item')

def testMappingSlicing():
	"""
	Test various slicing operations
	"""
	dest_dir = os.environ['PYTHONSCRATCH']
	source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

	mapping_base = HIVDatabase.MappingBase(source_dir, dest_dir,
									   'test_self')
	POSSIBLE_KEY = 'AJ302647.1'
	with open(seq_file) as handle:
		seq_val = SeqIO.parse(handle, 'fasta').next()
	test_key = POSSIBLE_KEY + seq_val.id

	mapping_base.AddtoShelf([seq_val])
	
	nose.tools.assert_true(mapping_base[test_key] != None,
							'Could not [] into MappingBase')
	test_mapping = mapping_base[test_key][0:5]
	nose.tools.assert_true(test_mapping != None,
							'Could not [:] into .is_match')

def testMappingBaseIter():
	"""
	Test the Iteration over MappingRecords
	"""
	
	dest_dir = os.environ['PYTHONSCRATCH']
	source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

	mapping_base = HIVDatabase.MappingBase(source_dir, dest_dir,
									   'test_self')
	POSSIBLE_KEY = 'AJ302647.1'
	with open(seq_file) as handle:
		seq_val = SeqIO.parse(handle, 'fasta').next()
	test_key = POSSIBLE_KEY + seq_val.id

	mapping_base.AddtoShelf([seq_val])

	for this_map in mapping_base.GetIter():
		nose.tools.assert_true(this_map != None, 
				'Could not generate individual mappings')
	
	for this_map in mapping_base.GetIter(WANTED_SUBTYPE = 'B'):
		nose.tools.assert_true(this_map != None, 
				'Could not generate individual mappings from a single subtype')
	
def testCheckSeqs():
	"""
	Test the Sequence Homology
	"""
	dest_dir = os.environ['PYTHONSCRATCH']
	source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

	mapping_base = HIVDatabase.MappingBase(source_dir, dest_dir,
									   'test_self')
	
	correct_seq = 'A'
	wrong_seq = 'QQQQQQ'
	
	nose.tools.assert_true(mapping_base.CheckSeqs(correct_seq) == 1,
							'Could not find "A" in all sequences.')
	nose.tools.assert_true(mapping_base.CheckSeqs(wrong_seq) == 0,
							'Found wrong_seq in a sequence')
	
def testCheckWindows():
	"""
	Test the FindWindows
	"""
	
	
	
def tearDownModule():
	time.sleep(5)
	file_list = os.listdir(os.environ['PYTHONSCRATCH'])
	for this_file in file_list:
		if this_file[-3:] != 'slf':
			os.remove(os.environ['PYTHONSCRATCH'] + this_file)

