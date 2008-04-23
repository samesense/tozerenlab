import os
from Bio import SeqIO
from numpy import *
import re
from pylab import *

os.chdir(os.curdir)

from PyELM import *


bkgHandle = open('test_data.fasta','rU')
bkgSeqs=[]

for thisSeq in SeqIO.parse(bkgHandle,'fasta'):
    bkgSeqs.append(thisSeq.seq.tostring())

bkgHandle.close()

ELMAnal = PyELM()

ELMAnal.ELMParser()
ELMAnal.LoadSeqs(bkgSeqs)

ELMAnal.ELMMatchSeqs()

counter=-1
for thisSeq in ELMAnal.GetSeqIterator():
    counter += 1
    print counter
    print sum(thisSeq)


