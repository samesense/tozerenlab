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

def MakeMultiCall(ELMData,SEQ_IND,WANTED_ELM,allBins):
    thisCall = []
    for i in xrange(len(allBins)-1):
        if sum(ELMData.ELMPosGen(SEQ_IND,allBins[i],allBins[i+1],WANTED_ELM)) > 0:
            thisCall.append(1)
        else:
            thisCall.append(0)
    return thisCall
    

def EvaluateBin(ELMData,WANTED_ELM,binEdges):
    thisBinVal=[0,0,0]
    for thisSeq in ELMAnal.GetSeqIterator(WANTED_ELM,(binEdges[0],binEdges[1])):
        CurrentBin=sum(thisSeq)
        if CurrentBin == 0:
            thisBinVal[0]+=1
        elif CurrentBin == 1:
            thisBinVal[1]+=1
        else:
            thisBinVal[2]+=1
    return (thisBinVal[0], thisBinVal[1], thisBinVal[2] )


def EvaluateAllBins(ELMData,WANTED_ELM,allBins):
    if not( tuple(allBins) in allBinDict ):

        TotalEmpty=0
        TotalCorrect=0
        TotalWrong=0
        for i in xrange(len(allBins)-1):
            if (allBins[i],allBins[i+1]) not in singlebinDict:
                singlebinDict[(allBins[i],allBins[i+1])] = EvaluateBin(ELMData,WANTED_ELM,(allBins[i],allBins[i+1]))
            
            TotalEmpty += singlebinDict[(allBins[i],allBins[i+1])][0]
            TotalCorrect += singlebinDict[(allBins[i],allBins[i+1])][1]
            TotalWrong += singlebinDict[(allBins[i],allBins[i+1])][2]
        funVal = (.25*TotalEmpty + TotalWrong) / (TotalCorrect + 1)
        allBinDict[tuple(allBins)]=(TotalEmpty,TotalCorrect,TotalWrong,funVal)

    return allBinDict[tuple(allBins)][3]


def MakeBreadth(globalQueue,allBins,BreadthCounter,DepthCounter):
    print (BreadthCounter,DepthCounter,EvaluateAllBins(ELMAnal,wanted_elm,allBins),globalQueue.qsize(),wanted_elm)
    if DepthCounter > MaxDepth:
        return

    if (BreadthCounter > MaxBreadth) & (allBins in DepthChecked):
        return
        
    LocalQueue = []
    LocalVals = []
        
    for k in xrange(len(allBins)-1):
        thisBinSet = deepcopy(allBins)
        newBinEdge=thisBinSet[k]+int((thisBinSet[k+1]-thisBinSet[k])/2)
        thisBinSet.insert(k+1,newBinEdge)
        LocalQueue.append(thisBinSet)
        LocalVals.append(EvaluateAllBins(ELMAnal,wanted_elm,thisBinSet))

    for k in xrange(1,len(allBins)-1):
        thisBinSet = deepcopy(allBins)
        thisBinSet.pop(k)
        LocalQueue.append(thisBinSet)
        LocalVals.append(EvaluateAllBins(ELMAnal,wanted_elm,thisBinSet))

    for k in xrange(1,len(allBins)-1):
        thisBinSet = deepcopy(allBins)
        newBinMove = int((thisBinSet[k+1]-thisBinSet[k])/2)
        LocalQueue.append(thisBinSet)
        LocalVals.append(EvaluateAllBins(ELMAnal,wanted_elm,thisBinSet))

    for k in xrange(1,len(allBins)-1):
        thisBinSet = deepcopy(allBins)
        newBinMove=int((thisBinSet[k]-thisBinSet[k-1])/2)
        thisBinSet[k]=thisBinSet[k]-newBinMove
        LocalQueue.append(thisBinSet)
        LocalVals.append(EvaluateAllBins(ELMAnal,wanted_elm,thisBinSet))

    BestInds=argsort(LocalVals)

    if BreadthCounter < MaxBreadth:
        for k in xrange(min(MaxWidth,len(BestInds))-1,0,-1):
            globalQueue.put(((LocalQueue[BestInds[k]],BreadthCounter+1,DepthCounter),LocalVals[BestInds[k]]))
    else:
        if (LocalVals[BestInds[0]] < EvaluateAllBins(ELMAnal,wanted_elm,allBins)):
            DepthChecked.append(allBins)
            globalQueue.put(((LocalQueue.pop(BestInds[0]),BreadthCounter+1,DepthCounter+1),LocalVals[BestInds[0]]))


    


def nonQueueWorker():
    GlobalCounter = 0
    while True:
        GlobalCounter += 1
        print GlobalCounter
        try:
            thisItem = BreadthQueue.get(True, 60)
            BreadthQueue.task_done()
        except:
            print 'queue empty'
            break
        
        if GlobalCounter < MaxIter:
            try:
                MakeBreadth(BreadthQueue,thisItem[0],thisItem[1],thisItem[2])
            except:
                print 'weird error'

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
