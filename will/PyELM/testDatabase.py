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

    mapping_base = HIVDatabase.MappingBase(source_dir, dest_dir,
                                           'test_self', 5)

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

	mapping_base = HIVDatabase.MappingBase(source_dir, dest_dir, 
											'test_self_KEEP')
	POSSIBLE_KEY = 'AJ302647.12_bkg_data_1'
	handle = open(seq_file)
	seq_iter = SeqIO.parse(handle, 'fasta')
	mapping_base.AddtoShelf(seq_iter)
	handle.close()

	nose.tools.assert_true(mapping_base.my_map_shelf.has_key(POSSIBLE_KEY),
							'Shelf is missing an Item')
	for this_key in mapping_base.my_map_shelf:
		match_sum = numpy.sum(mapping_base.my_map_shelf[this_key].is_match)
		f_look_sum = numpy.sum(mapping_base.my_map_shelf[this_key].f_look)
		nose.tools.assert_true(match_sum > 0, 'There were no matchings found')
		nose.tools.assert_true(f_look_sum > 0, 'F_look was not properly calculated')
	

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

	for this_key in mapping_base.my_map_shelf:
		match_sum = numpy.sum(mapping_base.my_map_shelf[this_key].is_match)
		f_look_sum = numpy.sum(mapping_base.my_map_shelf[this_key].f_look)
		nose.tools.assert_true(match_sum > 0, 'There were no matchings found')
		nose.tools.assert_true(f_look_sum > 0, 'F_look was not properly calculated')

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

							
def testAnnotGlobalHom():
	"""
	Test global homology figures
	"""
	dest_dir = os.environ['PYTHONSCRATCH']
	source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

	mapping_base = HIVDatabase.MappingBase(source_dir, dest_dir,
									   'test_self_KEEP')
									   
	for this_ref in mapping_base.ref_base:
		yield CheckAnnotGlobalHom, mapping_base, this_ref
		
def CheckAnnotGlobalHom(MAPPING_BASE, THIS_REF):
	
	dest_dir = os.environ['PYTHONSCRATCH'] + THIS_REF.seq_name + '_KEEP.pdf'
	MAPPING_BASE.WindowedHomology(THIS_REF.seq_name)
	
	len(THIS_REF.global_hom)
	len(THIS_REF.my_sequence)
	#nose.tools.assert_true(len(THIS_REF.global_hom) == len(THIS_REF.my_sequence),
	#						'Did not make Homology properly')
							
	THIS_REF.WriteGenesToDiagram()
	THIS_REF.AnnotGlobalHom()
	THIS_REF.this_genome.write(dest_dir, 'PDF')



# def testCheckWindows_SLOW():
	# """
	# Test the FindWindows
	# """
	# dest_dir = os.environ['PYTHONSCRATCH']
	# source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	# seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

	# mapping_base = HIVDatabase.MappingBase(source_dir, dest_dir,
									   # 'test_self_KEEP')
	# with open(dest_dir + 'output.slf', mode = 'w') as handle:
		# mapping_base.FindWindows(10, .4, handle)
	

def tearDownModule():
	checker = re.compile('.*KEEP.*')
	time.sleep(1)
	file_list = os.listdir(os.environ['PYTHONSCRATCH'])
	for this_file in file_list:
		if checker.match(this_file) == None:
			os.remove(os.environ['PYTHONSCRATCH'] + this_file)

    