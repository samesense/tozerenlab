import os
from Bio import SeqIO
from numpy import *
import re
from copy import *
from Queue import *
import cProfile
from threading import *
import pickle


os.chdir(os.curdir)

from PyELM import *
from pyQueueUtils import *

class SeqAns:
    def __init__(self):
        self.seq = []
        self.SimpleAns = []
        self.MultiAns = []
        self.ELM = []



bkgHandle = open('MedSeqs.pkl','rU')
bkgSeqs=pickle.load(bkgHandle)
##
##for thisSeq in SeqIO.parse(bkgHandle,'fasta'):
##    bkgSeqs.append(thisSeq.seq.tostring())
##
bkgHandle.close()




ELMAnal = PyELM()

##ELMAnal.ELMParser()
ELMAnal.AddELM('TESTCASE','CCCC')
ELMAnal.LoadSeqs(bkgSeqs.seq)

ELMAnal.ELMMatchSeqs()

ELMAnal.CalculateBins('TESTCASE')

##theseBins = [0, int(ELMAnal.MaxSize/2), ELMAnal.MaxSize]
##
##
##MaxBreadth = 4
##MaxDepth = 100
##MaxWidth = 5
##MaxIter = 1000
##
##
##
##
##for wanted_elm in ELMAnal.GetELMIterator():
##    BreadthQueue = PriorityQueue(-1)
##    print wanted_elm
##    BreadthQueue.put(((theseBins,0,0),0))
##       
##    
##    singlebinDict = {}
##    allBinDict = {}
##    DepthChecked = []
##    GlobalCounter = 0
##
##    nonQueueWorker()
##
####    for i in range(10):
####        t = Thread(target = nonQueueWorker)
####        t.setDaemon(True)
####        t.start()
####
####    BreadthQueue.join()
####
####
####
####
##    currentMin = 50000
##    currentMinInd = 0
##    for i in allBinDict:
##        if allBinDict[i][3] < currentMin:
##            currentMinInd = i
##            currentMin = allBinDict[currentMinInd][3]
##
##    print (wanted_elm, allBinDict[currentMinInd])
##
##    for i in xrange(len(ELMAnal.Sequences)-1):
##        print MakeMultiCall(ELMAnal,i,wanted_elm,currentMinInd)
##
##    
####
