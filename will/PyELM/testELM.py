import PyELM
import pickle
import os
import nose.tools
import sys
import re
import numpy

class SeqAns:
    """
    SeqAns
        A simple setup for keeping data and together for testing.
    """
    def __init__(self):
        self.seq = []
        self.SimpleAns = []
        self.MultiAns = []
        self.ELM = []

def AlignObs(shorterBins,longerBins):
    """
    AlignObs
        A function for aligning observations which are of different length
    """
    if len(shorterBins)==len(longerBins):
        matchVal=0
        for i in xrange(len(shorterBins)):
            if shorterBins[i]==shorterBins[i]:
                matchVal+=1
        return float(matchVal)/len(longerBins)
    
    delay = 0
    delayMax = len(longerBins) - len(shorterBins)
    pointerMat = []
    matchVal = 0
    for i in xrange(len(shorterBins)):
        if delay <= delayMax:
            if shorterBins[i] == longerBins[i+delay]:
                matchVal += 1
            elif shorterBins[i] == longerBins[i+delay+1]:
                matchVal +=1
                delay += 1      
    return float(matchVal)/len(longerBins)


def FileSetup(FILENAME):
    """
    FileSetup
        Loads and imports the SimpleAns data.
    """
    imported_data = SeqAns()
    file_loc = "C:\\Documents and Settings\\William Dampier\\My Documents\\PyELM\\"
    
    easy_handle = open(file_loc + FILENAME, 'r')
    data = pickle.load(easy_handle)
    easy_handle.close()

    imported_data.seq = data[0]
    imported_data.SimpleAns = data[1]
    imported_data.MultiAns = data[2]
    imported_data.ELM = data[3]

    return imported_data


def testSeqObject():
    """
    Test the SeqObject interface
    """

    test_dict = {}
    test_dict['TEST'] = re.compile('T')
    seq_obj = PyELM.SeqObject('AAAAAAAATAAA', ELM_DICT = test_dict)
    nose.tools.assert_not_equal(seq_obj, None)

def testGetResult():
    """
    Test window SeqObject searching
    """

    test_dict = {}
    test_dict['TEST'] = re.compile('T')
    seq_obj = PyELM.SeqObject('AAAAAAAATAAA', ELM_DICT = test_dict)

    nose.tools.assert_equal(seq_obj.GetResult('TEST', None), 1,
                            'Incorrect Count (No RANGE)')
    
    nose.tools.assert_equal(seq_obj.GetResult('TEST', (0, 6)), 0,
                            'Incorrect Count (with RANGE)')

    output = seq_obj.GetResult('TEST', None, OUTPUT_TYPE = 'array')


def testImport():
    """
    Test PyELM import
    """
    py_elm = PyELM.PyELM()
    nose.tools.assert_not_equal(py_elm, None)


def testImportSeq():
    """
    Test loading Sequence
    """
    file_list = ['EasyData.pkl', 'MedData.pkl', 'HardData.pkl']
    for this_file in file_list:
        yield CheckImportSeq, this_file

def CheckImportSeq(THIS_FILE):
    imported_data = FileSetup(THIS_FILE)
    py_elm = PyELM.PyELM()
    py_elm.LoadSeqs(imported_data.seq)


def testImportELM():
    """
    Test ELM importing
    """
    file_list = ['EasyData.pkl', 'MedData.pkl', 'HardData.pkl']
    for this_file in file_list:
        yield CheckImportELM, this_file

def CheckImportELM(THIS_FILE):
    imported_data = FileSetup(THIS_FILE)
    py_elm = PyELM.PyELM()
    py_elm.AddELM('TESTCASE', imported_data.ELM)
    py_elm.LoadSeqs(imported_data.seq)

def testParser():
    """
    Test ELM file parsing
    """
    py_elm = PyELM.PyELM()
    py_elm.ELMParser()

    nose.tools.assert_not_equal(py_elm.elm_dict, None)
    nose.tools.assert_not_equal(py_elm.elm_order, None)
    nose.tools.assert_not_equal(py_elm.elm_re_dict, None)


def CheckSimpleELM(THIS_FILE):

    imported_data = FileSetup(THIS_FILE)
    num_seqs = len(imported_data.seq)

    py_elm = PyELM.PyELM()
    py_elm.LoadSeqs(imported_data.seq)
    py_elm.AddELM('TESTCASE', imported_data.ELM)

    all_correct = 0    
    for i in xrange(num_seqs - 1):
        if py_elm.SimpleELMCall(i,'TESTCASE') == imported_data.SimpleAns[i]:
            all_correct += 1


    spec = float(all_correct) / num_seqs

    nose.tools.assert_true(spec > 0.9, msg = 'SimpleELM Spec: ' + str(spec))

def testSimpleELM():
    """
    Test SimpleELM call interface
    """
    file_list = ['EasyData.pkl', 'MedData.pkl', 'HardData.pkl']
    for this_file in file_list:
        yield CheckSimpleELM, this_file

def testIterator1():
    """
    Test iterator type 1
    """
    ITER_TYPE = 'ELM'
    ITER_RETURN = 'Name'
    ITER_RANGE = None

    elm_names = ['TESTCASE','TEST1','TEST2']
    py_elm = PyELM.PyELM()
    py_elm.AddELM(elm_names[0], 'TESTDATA')
    py_elm.AddELM(elm_names[1], 'BBBBB')
    py_elm.AddELM(elm_names[2], 'XXXXX')

    counter = 0
    for this_name in py_elm.GetIter(ITER_TYPE, ITER_RETURN, ITER_RANGE):
        nose.tools.assert_equal(this_name, elm_names[counter],
                                'ELM Name Iterator Failed')
        counter += 1

def testIterator2():
    """
    Test iterator type 2
    """
    ITER_TYPE = 'ELM'
    ITER_RETURN = 0
    ITER_RANGE = None

    imported_data = FileSetup('MedData.pkl')

    py_elm = PyELM.PyELM()
    py_elm.LoadSeqs(imported_data.seq)
    py_elm.AddELM('TESTCASE', imported_data.ELM)
    py_elm.AddELM('TEST1','BBBBB')
    py_elm.AddELM('TEST2','BBBBB')

    correct_ans = [True, False, False]
    counter = 0
    for this_val in py_elm.GetIter(ITER_TYPE, ITER_RETURN, ITER_RANGE):
        print this_val
        nose.tools.assert_equal(this_val,correct_ans[counter],
                                'ELM P/A Calls Failed')
        counter += 1
        

def testIterator3():
    """
    Test iterator type 3
    """
    ITER_TYPE = 'SEQ'
    ITER_RETURN = 'Sequence'

    imported_data = FileSetup('MedData.pkl')

    py_elm = PyELM.PyELM()
    py_elm.LoadSeqs(imported_data.seq)

    counter = 0
    for this_val in py_elm.GetIter(ITER_TYPE, ITER_RETURN, None):
        nose.tools.assert_equal(len(this_val),
                                len(imported_data.seq[counter]),
                                'Sequence Iteration Failed')
        counter += 1

    for this_val in py_elm.GetIter(ITER_TYPE, ITER_RETURN, [0, 200]):
        nose.tools.assert_equal(len(this_val),
                                200, 'Sequence Iteration Failed')
        
def testIterator4():
    """
    Test iterator type 4
    """
    ITER_TYPE = 'SEQ'
    ITER_RETURN = 'TESTCASE'
    ITER_RANGE = None
    
    imported_data = FileSetup('MedData.pkl')

    py_elm = PyELM.PyELM()
    py_elm.AddELM('TESTCASE', imported_data.ELM)
    py_elm.LoadSeqs(imported_data.seq)


    counter = 0
    for this_val in py_elm.GetIter(ITER_TYPE, ITER_RETURN, ITER_RANGE):
        nose.tools.assert_equal(this_val, imported_data.SimpleAns[counter] == 1,
                                'Sequence Iteration Failed')
        counter += 1


def testMultiELM():
    """
    Test MultiELM
    """
    file_list = ['EasyData.pkl', 'MedData.pkl', 'HardData.pkl']
    for this_file in file_list:
        yield CheckMultiELM, this_file

def CheckMultiELM(THIS_FILE):

    imported_data = FileSetup(THIS_FILE)
    num_seqs = len(imported_data.seq)

    py_elm = PyELM.PyELM()
    py_elm.force_bins = []
    py_elm.LoadSeqs(imported_data.seq)
    py_elm.AddELM('TESTCASE', imported_data.ELM)

    total_val = float(0)
    py_elm.CalculateBins('TESTCASE')
    for i in xrange(num_seqs - 1):
        this_ans = py_elm.MultiELMCall(i, 'TESTCASE')
        if len(imported_data.MultiAns[i]) <= len(this_ans):
            total_val += AlignObs(imported_data.MultiAns[i], this_ans)
        else:
            total_val += AlignObs(this_ans, imported_data.MultiAns[i])

    spec = total_val / num_seqs

    nose.tools.assert_true(spec > 0.8, msg = 'SimpleELM Spec: ' + str(spec))