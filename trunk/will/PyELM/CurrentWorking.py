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



##bkgHandle = open('bkg_aa_seq.fasta','rU')
##bkgSeqs=[]
####bkgSeqs=pickle.load(bkgHandle)
##
##for thisSeq in SeqIO.parse(bkgHandle,'fasta'):
##    bkgSeqs.append(thisSeq.seq.tostring())
##
##bkgHandle.close()


prevRun = open('backupData8.pkl','r')

ELMAnal = pickle.load(prevRun)

prevRun.close()

##ELMAnal = PyELM()
##
##ELMAnal.ELMParser()
##ELMAnal.LoadSeqs(bkgSeqs)
##
##ELMAnal.ELMMatchSeqs()
##ELMAnal.CalcTimeOut=1000
print 'Calculating Bins'

for wanted_elm in ELMAnal.GetELMIterator():
    print wanted_elm
    if wanted_elm in ELMAnal.PreCalcELMBinDict:
        print 'Already Calculated'
        continue

    try:
        ELMAnal.CalculateBins(wanted_elm)
    finally:
        print 'Backing Up'
        backupHandle = open('backupData9.pkl','w')
        pickle.dump(ELMAnal,backupHandle)
        backupHandle.close()
    
