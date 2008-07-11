from __future__ import with_statement
import nose.tools

from testELM import SeqAns, FileSetup
from Bio import SeqIO
import itertools
import os
import PyVirus

def testLoading():
    """
    Test Loading the Virus classes
    """
    possible_classes = [PyVirus.ViralSeq, PyVirus.BkgSeq,
                        PyVirus.PatSeq]
    file_loc = os.environ['MYDOCPATH'] + 'PyELM\\'
    with open(file_loc + '50_seqs.fasta', mode = 'r') as file_handle:
        this_iter = itertools.izip(SeqIO.parse(file_handle, 'fasta'),
                                   iter(possible_classes))
        for this_test in this_iter:
            yield CheckLoading, this_test[1], this_test[0]
        
def CheckLoading(INPUT_CLASS, INPUT_REC):
    genome = INPUT_CLASS(INPUT_REC.seq.tostring(), 'testFile')
    nose.tools.assert_not_equal(genome, None)
        

def testGenotyping():
    """
    Test Genotyping the Virus classes
    """
    possible_classes = [PyVirus.ViralSeq, PyVirus.BkgSeq,
                        PyVirus.PatSeq]
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
    base_dir = os.environ['MYDOCPATH'] + 'hivsnppred\\HIVRefs\\'
    for this_file in filter(lambda x: x[-3:] == '.gb', os.listdir(base_dir)):
        yield CheckRefLoading, base_dir+this_file
    

def CheckRefLoading(INPUT_FILE):
    genome = PyVirus.RefSeq(INPUT_FILE)
    nose.tools.assert_not_equal(genome.annotation, None)

def testRefBaseLoading():
    """
    Test RefBase loading
    """
    base_dir = os.environ['MYDOCPATH'] + 'hivsnppred\\HIVRefs\\'
    dest_dir = 'C:\\local_blast\\PyELMData\\'
    ref_base = PyVirus.RefBase(base_dir, dest_dir)
    nose.tools.assert_not_equal(ref_base, None, 'Could not Load Data')
    nose.tools.assert_not_equal(ref_base.ref_seqs, None, 'Data not present')
