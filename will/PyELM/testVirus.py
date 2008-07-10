from __future__ import with_statement
import nose.tools
import PyVirus
from testELM import SeqAns, FileSetup
from Bio import SeqIO
import itertools
import os


def testLoading():
    """
    Test Loading the Virus classes
    """
    possible_classes = [PyVirus.ViralSeq, PyVirus.BkgSeq,
                        PyVirus.PatSeq]
    file_loc = "C:\\Documents and Settings\\Will\\My Documents\\PyELM\\"
    with open(file_loc + '50_seqs.fasta', mode = 'r') as file_handle:
        this_iter = itertools.izip(SeqIO.parse(file_handle, 'fasta'),
                                   iter(possible_classes))
        for this_test in this_iter:
            yield CheckLoading, this_test[1], this_test[0]
        
def CheckLoading(INPUT_CLASS, INPUT_REC):
    genome = INPUT_CLASS(INPUT_REC.seq.tostring())
    nose.tools.assert_not_equal(genome, None)
        

def testGenotyping():
    """
    Test Genotyping the Virus classes
    """
    possible_classes = [PyVirus.ViralSeq, PyVirus.BkgSeq,
                        PyVirus.PatSeq]
    file_loc = "C:\\Documents and Settings\\Will\\My Documents\\PyELM\\"
    with open(file_loc + '50_seqs.fasta', mode = 'r') as file_handle:
        this_iter = itertools.izip(SeqIO.parse(file_handle, 'fasta'),
                                   iter(possible_classes))
        for this_test in this_iter:
            yield CheckGenotyping, this_test[1], this_test[0]
        
def CheckGenotyping(INPUT_CLASS, INPUT_REC):
    genome = INPUT_CLASS(INPUT_REC.seq.tostring())
    genome.DetSubtype()
    nose.tools.assert_not_equal(genome, None)


def testRefLoading():
    """
    Test Genbank parsing of Reference Sequences
    """
    base_dir = 'C:\\Documents and Settings\\Will\\My Documents\\hivsnppredsvn\\HIVRefs\\'
    for this_file in filter(lambda x: x[-3:] == '.gb', os.listdir(base_dir)):
        yield CheckRefLoading, base_dir+this_file
    

def CheckRefLoading(INPUT_FILE):
    genome = PyVirus.RefSeq(INPUT_FILE)
    nose.tools.assert_not_equal(genome.annotation, None)
