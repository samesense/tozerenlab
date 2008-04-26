import unittest
from PyELM import *
import pickle
import os

os.chdir(os.curdir)


class SeqAns:
    def __init__(self):
        self.seq = []
        self.SimpleAns = []
        self.MultiAns = []
        self.ELM = []


class TestEasyCase(unittest.TestCase):
    def setUp(self):
        try:
            EasyHandle = open('EasySeqs.pkl','r')
            self.EasyData = pickle.load(EasyHandle)
            self.Anal = PyELM()
            self.Anal.LoadSeqs(self.EasyData.seq)
            self.Anal.AddELM('TESTCASE',self.EasyData.ELM)
            
        finally:
            EasyHandle.close()

    def tearDown(self):
        self.EasyData = None

    def testSimpleELM(self):
        AllCorrect = True
        for i in xrange(len(self.EasyData.seq)-1):
            if self.Anal.SimpleELMCall(i,'TESTCASE') != self.EasyData.SimpleAns[i]:
                AllCorrect = False

        self.failUnless(AllCorrect,'Could not Predict SimpleELM.')

    def testMultiELM(self):
        AllCorrect = True
        for i in xrange(len(self.EasyData.seq)-1):
            if self.Anal.MultiELMCall(i,'TESTCASE') != self.EasyData.MultiAns[i]:
                AllCorrect = False

        self.failUnless(AllCorrect,'Could not Predict SimpleELM.')
