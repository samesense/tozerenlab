import os
from Bio import SeqIO
from numpy import *
import re
from copy import *
from Queue import *

HOME_DIR="C:\\Documents and Settings\\William Dampier\\My Documents\\ELM_Motif_finder\\ELM_RAW_DOWNLOAD\\"

os.chdir(os.curdir)

from PyELM import *

TestedSolutions = {}


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
    if not( tuple(allBins) in TestedSolutions ):
        TotalEmpty=0
        TotalCorrect=0
        TotalWrong=0
        for i in xrange(len(allBins)-1):
            if (allBins[i],allBins[i+1]) not in CalculatedBins:
                CalculatedBins[(allBins[i],allBins[i+1])] = EvaluateBin(ELMData,WANTED_ELM,(allBins[i],allBins[i+1]))
            
            TotalEmpty += CalculatedBins[(allBins[i],allBins[i+1])][0]
            TotalCorrect += CalculatedBins[(allBins[i],allBins[i+1])][1]
            TotalWrong += CalculatedBins[(allBins[i],allBins[i+1])][2]
        TestedSolutions[tuple(theseBins)]=(TotalEmpty,TotalCorrect,TotalWrong)
    return TestedSolutions[tuple(theseBins)][2]

bkgHandle = open('test_data.fasta','rU')
bkgSeqs=[]

for thisSeq in SeqIO.parse(bkgHandle,'fasta'):
    bkgSeqs.append(thisSeq.seq.tostring())

bkgHandle.close()

ELMAnal = PyELM()

ELMAnal.ELMParser(HOME_DIR)
ELMAnal.LoadSeqs(bkgSeqs)

ELMAnal.ELMMatchSeqs()

theseBins = [0, int(ELMAnal.MaxSize/2), ELMAnal.MaxSize]

wanted_elm='CLV_PCSK_KEX2_1'

CalculatedBins = {}


##TotalEmpty=0
##TotalCorrect=0
##TotalWrong=0
##for i in xrange(len(theseBins)-1):
##    if (theseBins[i],theseBins[i+1]) not in CalculatedBins:
##        CalculatedBins[(theseBins[i],theseBins[i+1])] = EvaluateBin(ELMAnal,wanted_elm,(theseBins[i],theseBins[i+1]))
##    TotalEmpty += CalculatedBins[(theseBins[i],theseBins[i+1])][0]
##    TotalCorrect += CalculatedBins[(theseBins[i],theseBins[i+1])][1]
##    TotalWrong += CalculatedBins[(theseBins[i],theseBins[i+1])][2]
##TestedSolutions[tuple(theseBins)]=(TotalEmpty,TotalCorrect,TotalWrong)
##
##CurrentBest=TestedSolutions[tuple(theseBins)][2]
##BestMethod=''
##BestInd=0

##for p in xrange(40):
##    #try adding Edges in the middle of current edges
##    for k in xrange(len(theseBins)-1):
##
##        newBinEdge=theseBins[k]+int((theseBins[k+1]-theseBins[k])/2)
##        theseBins.insert(k+1,newBinEdge)
##        if tuple(theseBins) in TestedSolutions:
##            continue
##
##        EvaluateAllBins(ELMAnal,wanted_elm,theseBins)
##
##        if TestedSolutions[tuple(theseBins)][2] < CurrentBest:
##            BestMethod='ADD'
##            BestInd=k
##            CurrentBest=TestedSolutions[tuple(theseBins)][2]
##
##        theseBins.pop(k+1)
##
##    #try moving edges to the midpoint
##    for k in xrange(1,len(theseBins)-1):
##
##        newBinMove=int((theseBins[k+1]-theseBins[k])/2)
##        theseBins[k]=theseBins[k]+newBinMove
##
##        if tuple(theseBins) in TestedSolutions:
##            continue        
##
##        EvaluateAllBins(ELMAnal,wanted_elm,theseBins)
##
##        if TestedSolutions[tuple(theseBins)][2] < CurrentBest:
##            BestMethod='MoveRight'
##            BestInd=k
##            CurrentBest=TestedSolutions[tuple(theseBins)][2]
##
##        theseBins[k]=theseBins[k]-newBinMove
##
##    for k in xrange(1,len(theseBins)-1):
##        newBinMove=int((theseBins[k]-theseBins[k-1])/2)
##        
##        theseBins[k]=theseBins[k]-newBinMove
##
##        
##        if tuple(theseBins) in TestedSolutions:
##            continue        
##
##        EvaluateAllBins(ELMAnal,wanted_elm,theseBins)
##
##        if TestedSolutions[tuple(theseBins)][2] < CurrentBest:
##            BestMethod='MoveLeft'
##            BestInd=k
##            CurrentBest=TestedSolutions[tuple(theseBins)][2]
##
##        theseBins[k]=theseBins[k]+newBinMove
##
##
##
##    if BestMethod=='ADD':
##        newBinEdge=theseBins[BestInd]+int((theseBins[BestInd+1]-theseBins[BestInd])/2)
##        theseBins.insert(BestInd+1,newBinEdge)
##    elif BestMethod == 'MoveRight':
##        newBinMove=int((theseBins[k+1]-theseBins[k])/2)
##        theseBins[k]=theseBins[k]+newBinMove
##    elif BestMethod == 'MoveLeft':
##        newBinMove=int((theseBins[k]-theseBins[k-1])/2)
##        theseBins[k]=theseBins[k]-newBinMove
##


globalQueue = Queue(-1)
MaxBreadth = 3
MaxDepth = 5

def MakeBreadth(allBins,BreadthCounter,DepthCounter):
    
    print (BreadthCounter,DepthCounter,globalQueue.qsize())

    if DepthCounter > MaxDepth:
        return
    
    if BreadthCounter > MaxBreadth:
        print 'CurrentVals'
        LocalQueue = []
        
    for k in xrange(len(allBins)-1):
        thisBinSet = deepcopy(allBins)
        newBinEdge=thisBinSet[k]+int((thisBinSet[k+1]-thisBinSet[k])/2)
        thisBinSet.insert(k+1,newBinEdge)
        if BreadthCounter > MaxBreadth:
            LocalQueue.append(thisBinSet)
        else:
            globalQueue.put((thisBinSet,BreadthCounter+1,DepthCounter))

    for k in xrange(1,len(allBins)-1):
        thisBinSet = deepcopy(allBins)
        thisBinSet.pop(k)
        if BreadthCounter > MaxBreadth:
            LocalQueue.append(thisBinSet)
        else:
            globalQueue.put((thisBinSet,BreadthCounter+1,DepthCounter))

    for k in xrange(1,len(allBins)-1):
        thisBinSet = deepcopy(allBins)
        newBinMove = int((thisBinSet[k+1]-thisBinSet[k])/2)
        if BreadthCounter > MaxBreadth:
            LocalQueue.append(thisBinSet)
        else:
            globalQueue.put((thisBinSet,BreadthCounter+1,DepthCounter))

    for k in xrange(1,len(allBins)-1):
        thisBinSet = deepcopy(allBins)
        newBinMove=int((thisBinSet[k]-thisBinSet[k-1])/2)
        thisBinSet[k]=thisBinSet[k]-newBinMove
        if BreadthCounter > MaxBreadth:
            LocalQueue.append(thisBinSet)
        else:
            globalQueue.put((thisBinSet,BreadthCounter+1,DepthCounter))

    if BreadthCounter > MaxBreadth:
        bestBinSet = LocalQueue.pop()
        bestVal = EvaluateAllBins(ELMAnal,wanted_elm,bestBinSet)
        while len(LocalQueue) !=0:
            thisBin = LocalQueue.pop()
            thisVal = EvaluateAllBins(ELMAnal,wanted_elm,thisBin)
            if thisVal < bestVal:
                bestVal = thisVal
                bestBinSet = deepcopy(thisBin)

        globalQueue.put((bestBinSet,BreadthCounter,DepthCounter+1))
            

globalQueue.put((theseBins,0,0))
while not(globalQueue.empty()):
    thisItem = globalQueue.get()
    globalQueue.task_done()
    MakeBreadth(thisItem[0],thisItem[1],thisItem[2])










