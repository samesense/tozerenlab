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
    seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

    mapping_base = HIVDatabase.MappingBase(source_dir, seq_file,
                                           dest_dir, 'test_self.slf')

    nose.tools.assert_not_equal(mapping_base, None)
    

def testMappingBase():
    """
    Test loading MappingBase
    """
    dest_dir = os.environ['PYTHONSCRATCH']
    source_dir = os.environ['MYDOCPATH'] + 'hivsnppred\\HIVRefs\\'
    seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

    mapping_base = HIVDatabase.MappingBase(source_dir, seq_file,
                                           dest_dir, 'test_self.slf')

    nose.tools.assert_not_equal(mapping_base, None)
    
def testCreateShelf():
    """
    Test CreateShelf
    """
    dest_dir = os.environ['PYTHONSCRATCH']
    source_dir = os.environ['MYDOCPATH'] + 'hivsnppred\\HIVRefs\\'
    seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

    mapping_base = HIVDatabase.MappingBase(source_dir, seq_file,
                                           dest_dir, 'test_self.slf')

    mapping_base.CreateShelf()
    nose.tools.assert_not_equal(len(mapping_base.my_shelf.keys()), 0,
                                'Shelf has no items')



def tearDownModule():
    file_list = os.listdir(os.environ['PYTHONSCRATCH'])
    for this_file in file_list:
        os.remove(os.environ['PYTHONSCRATCH'] + this_file)

