from __future__ import with_statement
import nose.tools
import HIVDatabase
import PyVirus
import os
from Bio import SeqIO
import itertools

def setupModule():
    dest_dir = os.environ['PYTHONSCRATCH']
    source_dir = os.environ['MYDOCPATH'] + 'hivsnppred\\HIVRefs\\'
    
    PyVirus.RefBase(source_dir, dest_dir, BUILD = True)



def testMappingRecord():
    """
    Test loading MappingRecord
    """
    dest_dir = os.environ['PYTHONSCRATCH']
    source_dir = os.environ['MYDOCPATH'] + 'hivsnppred\\HIVRefs\\'
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
    source_dir = os.environ['MYDOCPATH'] + 'hivsnppred\\HIVRefs\\'

    mapping_base = HIVDatabase.MappingBase(source_dir, dest_dir,
                                           'test_self.slf')

    nose.tools.assert_not_equal(mapping_base, None)
    

def testMappingBase():
    """
    Test loading MappingBase
    """
    dest_dir = os.environ['PYTHONSCRATCH']
    source_dir = os.environ['MYDOCPATH'] + 'hivsnppred\\HIVRefs\\'
    seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

    mapping_base = HIVDatabase.MappingBase(source_dir, 
                                           dest_dir, 'test_self.slf')
    
    nose.tools.assert_not_equal(mapping_base, None)
    
def testAddtoShelf():
    """
    Test AddtoShelf
    """
    dest_dir = os.environ['PYTHONSCRATCH']
    source_dir = os.environ['MYDOCPATH'] + 'hivsnppred\\HIVRefs\\'
    seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

    mapping_base = HIVDatabase.MappingBase(source_dir, dest_dir,
                                           'test_self.slf')
    POSSIBLE_KEY = 'AJ302647.1'
    with open(seq_file) as handle:
        seq_iter = itertools.izip(SeqIO.parse(handle, 'fasta'),
                                  itertools.repeat(None, 5))

        for this_seq in seq_iter:
            mapping_base.AddtoShelf([this_seq[0]])
            last_id = this_seq[0].id
    test_key = POSSIBLE_KEY + last_id
    nose.tools.assert_true(mapping_base.my_shelf.has_key(test_key),
                                'Shelf is missing an Item')

def testSaveShelf():
    """
    Test SaveShelf and autoloading
    """

    dest_dir = os.environ['PYTHONSCRATCH']
    source_dir = os.environ['MYDOCPATH'] + 'hivsnppred\\HIVRefs\\'
    seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

    mapping_base = HIVDatabase.MappingBase(source_dir, dest_dir,
                                           'test_self.slf')
    POSSIBLE_KEY = 'AJ302647.1'
    with open(seq_file) as handle:
        seq_val = SeqIO.parse(handle, 'fasta').next()
    test_key = POSSIBLE_KEY + seq_val.id

    mapping_base.AddtoShelf([seq_val])
    mapping_base.SaveShelf()

    del(mapping_base)

    mapping_base = HIVDatabase.MappingBase(source_dir, dest_dir,
                                           'test_self.slf')
    
    nose.tools.assert_true(mapping_base.my_shelf.has_key(test_key),
                                'Shelf is missing an Item')


def tearDownModule():
    file_list = os.listdir(os.environ['PYTHONSCRATCH'])
    for this_file in file_list:
        os.remove(os.environ['PYTHONSCRATCH'] + this_file)

