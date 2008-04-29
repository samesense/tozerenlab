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

def AlignObs(shorterBins,longerBins):
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

def suite():
    suite = unittest.TestSuite()
    suite.addTest(TestEasyCase('testSimpleELM'))
    suite.addTest(TestEasyCase('testMultiELM'))
    suite.addTest(TestMedCase('testSimpleELM'))
    suite.addTest(TestMedCase('testMultiELM'))
    suite.addTest(TestHardCase('testSimpleELM'))
    suite.addTest(TestHardCase('testMultiELM'))
    return suite

def RunTest():
    thissuite=suite()
    unittest.TextTestRunner(verbosity=2).run(thissuite)

class TestEasyCase(unittest.TestCase):
    def setUp(self):
        try:
            EasyHandle = open('EasySeqs.pkl','r')
            self.EasyData = pickle.load(EasyHandle)
            self.Anal = PyELM()
            self.Anal.LoadSeqs(self.EasyData.seq)
            self.Anal.AddELM('TESTCASE',self.EasyData.ELM)
            self.Anal.CalcTimeOut=100
            
        finally:
            EasyHandle.close()

    def tearDown(self):
        self.EasyData = None

    def testSimpleELM(self):
        AllCorrect = 0
        for i in xrange(len(self.EasyData.seq)-1):
            if self.Anal.SimpleELMCall(i,'TESTCASE') == self.EasyData.SimpleAns[i]:
                AllCorrect += 1

        spec = float(AllCorrect)/len(self.Anal.Sequences)

        self.failIf(spec < 0.6, 'SinlgeELM Prediction failed at 0.6')
        self.failIf(spec < 0.7, 'SinlgeELM Prediction failed at 0.7')
        self.failIf(spec < 0.8, 'SinlgeELM Prediction failed at 0.8')
        self.failIf(spec < 0.9, 'SinlgeELM Prediction failed at 0.9')
        self.failIf(spec < 0.99, 'SinlgeELM Prediction failed at 0.99')

    def testMultiELM(self):
        totalVal = float(0)
        self.Anal.CalculateBins('TESTCASE')
        
        for i in xrange(len(self.EasyData.seq)-1):
            thisAns = self.Anal.MultiELMCall(i,'TESTCASE')
            if len(self.EasyData.MultiAns[i]) <= len(thisAns):
                totalVal += AlignObs(self.EasyData.MultiAns[i],thisAns)
            else:
                totalVal += AlignObs(thisAns,self.EasyData.MultiAns[i])

        spec = totalVal/len(self.EasyData.seq)

        self.failIf(spec < 0.6, 'MultiELM Prediction failed at 0.6')
        self.failIf(spec < 0.7, 'MultiELM Prediction failed at 0.7')
        self.failIf(spec < 0.8, 'MultiELM Prediction failed at 0.8')
        self.failIf(spec < 0.9, 'MultiELM Prediction failed at 0.9')
        self.failIf(spec < 0.99, 'MultiELM Prediction failed at 0.99')
        

   
class TestMedCase(unittest.TestCase):
    def setUp(self):
        try:
            MedHandle = open('MedSeqs.pkl','r')
            self.MedData = pickle.load(MedHandle)
            self.Anal = PyELM()
            self.Anal.LoadSeqs(self.MedData.seq)
            self.Anal.AddELM('TESTCASE',self.MedData.ELM)
            self.Anal.CalcTimeOut=100
        finally:
            MedHandle.close()

    def tearDown(self):
        self.MedData = None

    def testSimpleELM(self):
        AllCorrect = 0
        for i in xrange(len(self.MedData.seq)-1):
            if self.Anal.SimpleELMCall(i,'TESTCASE') == self.MedData.SimpleAns[i]:
                AllCorrect += 1

        spec = float(AllCorrect)/len(self.Anal.Sequences)
        
        self.failIf(spec < 0.6, 'SinlgeELM Prediction failed at 0.6')
        self.failIf(spec < 0.7, 'SinlgeELM Prediction failed at 0.7')
        self.failIf(spec < 0.8, 'SinlgeELM Prediction failed at 0.8')
        self.failIf(spec < 0.9, 'SinlgeELM Prediction failed at 0.9')
        self.failIf(spec < 0.99, 'SinlgeELM Prediction failed at 0.99')

    def testMultiELM(self):
        totalVal = float(0)
        self.Anal.CalculateBins('TESTCASE') 
        for i in xrange(len(self.MedData.seq)-1):
            thisAns = self.Anal.MultiELMCall(i,'TESTCASE')
            if len(self.MedData.MultiAns[i]) <= len(thisAns):
                totalVal += AlignObs(self.MedData.MultiAns[i],thisAns)
            else:
                totalVal += AlignObs(thisAns,self.MedData.MultiAns[i])

        spec = totalVal/len(self.MedData.seq)

        self.failIf(spec < 0.6, 'MultiELM Prediction failed at 0.6')
        self.failIf(spec < 0.7, 'MultiELM Prediction failed at 0.7')
        self.failIf(spec < 0.8, 'MultiELM Prediction failed at 0.8')
        self.failIf(spec < 0.9, 'MultiELM Prediction failed at 0.9')
        self.failIf(spec < 0.99, 'MultiELM Prediction failed at 0.99')
        

       

class TestHardCase(unittest.TestCase):
    def setUp(self):
        try:
            MedHandle = open('HardSeqs.pkl','r')
            self.MedData = pickle.load(MedHandle)
            self.Anal = PyELM()
            self.Anal.LoadSeqs(self.MedData.seq)
            self.Anal.AddELM('TESTCASE',self.MedData.ELM)
            self.Anal.CalcTimeOut=100
            
        finally:
            MedHandle.close()

    def tearDown(self):
        self.MedData = None

    def testSimpleELM(self):
        AllCorrect = 0
        for i in xrange(len(self.MedData.seq)-1):
            if self.Anal.SimpleELMCall(i,'TESTCASE') == self.MedData.SimpleAns[i]:
                AllCorrect += 1

        spec = float(AllCorrect)/len(self.Anal.Sequences)

        self.failIf(spec < 0.6, 'SinlgeELM Prediction failed at 0.6')
        self.failIf(spec < 0.7, 'SinlgeELM Prediction failed at 0.7')
        self.failIf(spec < 0.8, 'SinlgeELM Prediction failed at 0.8')
        self.failIf(spec < 0.9, 'SinlgeELM Prediction failed at 0.9')
        self.failIf(spec < 0.99, 'SinlgeELM Prediction failed at 0.99')


        self.failUnless(AllCorrect,'Could not Predict SimpleELM.')

    def testMultiELM(self):
        totalVal = float(0)
        self.Anal.CalculateBins('TESTCASE')
        for i in xrange(len(self.MedData.seq)-1):
            #print 'Guessed: ', self.Anal.MultiELMCall(i,'TESTCASE')
            #print 'Known: ', self.MedData.MultiAns[i]
            thisAns = self.Anal.MultiELMCall(i,'TESTCASE')
            if len(self.MedData.MultiAns[i]) <= len(thisAns):
                totalVal += AlignObs(self.MedData.MultiAns[i],thisAns)
            else:
                totalVal += AlignObs(thisAns,self.MedData.MultiAns[i])

        spec = totalVal/len(self.MedData.seq)

        self.failIf(spec < 0.6, 'MultiELM Prediction failed at 0.6')
        self.failIf(spec < 0.7, 'MultiELM Prediction failed at 0.7')
        self.failIf(spec < 0.8, 'MultiELM Prediction failed at 0.8')
        self.failIf(spec < 0.9, 'MultiELM Prediction failed at 0.9')
        self.failIf(spec < 0.99, 'MultiELM Prediction failed at 0.99')
        


