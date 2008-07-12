from __future__ import with_statement
import nose.tools
from testELM import SeqAns, FileSetup
from Bio import SeqIO
import itertools
import os
import PyVirus
import tempfile
import time


def testLoading():
    """
    Test Loading the Virus classes
    """
    possible_classes = [PyVirus.ViralSeq, PyVirus.BkgSeq,
                        PyVirus.PatSeq]
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
    base_dir = os.environ['MYDOCPATH'] + 'hivsnppred\\HIVRefs\\'
    dest_dir = os.environ['PYTHONSCRATCH']
    ref_base = PyVirus.RefBase(base_dir, dest_dir)
    nose.tools.assert_not_equal(ref_base, None, 'Could not Load Data')
    nose.tools.assert_not_equal(ref_base.ref_seqs, None, 'Data not present')


def testGenomeToBioPython():
    """
    Test the whole genome conversion to BioPython SeqRecord
    """
    possible_classes = [PyVirus.ViralSeq, PyVirus.BkgSeq,
                        PyVirus.PatSeq]
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

    correct_suffix_list = ['hr', 'in', 'sd', 'si', 'sq']
    correct_base_list = ['ref_aa.fasta.p', 'ref_nt.fasta.n']
    base_dir = os.environ['MYDOCPATH'] + 'hivsnppred\\HIVRefs\\'
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
    base_dir = os.environ['MYDOCPATH'] + 'hivsnppred\\HIVRefs\\'
    dest_dir = os.environ['PYTHONSCRATCH']
    ref_base = PyVirus.RefBase(base_dir, dest_dir)
    possible_classes = [PyVirus.ViralSeq, PyVirus.BkgSeq,
                        PyVirus.PatSeq]
    file_loc = os.environ['MYDOCPATH'] + 'PyELM\\'
    with open(file_loc + '50_seqs.fasta', mode = 'r') as file_handle:
        this_iter = itertools.izip(SeqIO.parse(file_handle, 'fasta'),
                                   iter(possible_classes))
        for this_test in this_iter:
            yield CheckTranslateAll, ref_base, this_test[1], this_test[0]

def CheckTranslateAll(REF_BASE, INPUT_CLASS, INPUT_RECORD):
    
    genome = INPUT_CLASS(INPUT_RECORD.seq.tostring(), 'testseq')
    genome.TranslateAll(REF_BASE)

    nose.tools.assert_true(genome.annotation.has_key('env'))

    


def tearDownModule():
    file_list = os.listdir(os.environ['PYTHONSCRATCH'])
    for this_file in file_list:
        os.remove(os.environ['PYTHONSCRATCH'] + this_file)



















    
