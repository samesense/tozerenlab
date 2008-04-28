import os
from numpy import *
import re
from pyQueueUtils import *
from copy import *


class PyELM:

    def __init__(self):

        self.ELMDict={}
        self.ELMOrder=[]
        self.Sequences=[]
        self.ELMreDict={}
        self.MaxSize=0

        self.__FoundCorrect = False

        self.PreCalchistArray=[]        
        self.PreCalcELMMatchSeqs=[]
        self.PreCalcELMBinDict={}

        self.__singlebinDict={}
        self.__allBinDict={}
        self.__DepthChecked=[]
        self.__WorkingQueue=PriorityQueue(-1)

        
        self.MaxBreadth = 4
        self.MaxDepth = 100
        self.MaxWidth = 5
        self.MaxIter = 1000

    def ELMParser(self,DIRECTORY="C:\Documents and Settings\Will\My Documents\ELM_Motif_finder\ELM_RAW_DOWNLOAD\\"):
        """
        Parser
            Parses a directory of HTML files to extract the ELM names and Regular Expressions and instanciates it into SELF.

            Parser(DIRECTORY)

            DIRECTORY       The directory containing the HTML files downloaded from the ELM database.
                            If left empty then C:\Documents and Settings\Will\My Documents\PyELM\ELM_RAW_DOWNLOAD\

        """


        fileList=os.listdir(DIRECTORY)
        allRegExps={}
        
        for i in xrange(1,len(fileList)-1):
            
            thisFile=open(DIRECTORY+fileList[i],'r')
            ELMname=fileList[i].rpartition('.')
            self.ELMOrder.append(ELMname[0])
            currentLine=''
            
            while currentLine.find('Pattern:') == -1:
                currentLine=thisFile.readline()
            
            currentLine=thisFile.readline()
            
            allRegExps[ELMname[0]]=currentLine[currentLine.find('>')+1:currentLine.find('/')-2]
            thisFile.close()
        
        self.ELMDict=allRegExps

    def AddELM(self,Name,RegExp):
        """ 
        AddELM
            Loads a single ELM RegExp into the database.

        AddELM(Name,RegExp)

        """
        self.ELMOrder.append(Name)
        self.ELMDict[Name] = RegExp
        self.ELMreDict[Name] = re.compile(RegExp)
        self.PreCalcELMMatchSeqs=[]
        self.PreCalchistArray=[]
        
    def LoadSeqs(self,SEQUENCES):
        """
        LoadSeqs
            Loads a set of sequences into this instance for further analysis.

            LoadSeqs(SEQUENCES)

            SEQUENCES       A LIST of characters representing Amino Acids


        """

        self.Sequences = SEQUENCES
        self.MaxSize=0
        for i in xrange(len(self.Sequences)):
            if len(self.Sequences[i]) > self.MaxSize :
                self.MaxSize=len(self.Sequences[i])


    

    def ELMMatchSeqs(self):
        """
        ELMMatchSeqs
            Searches for the regular expressions of ELM_DICT within the SEQUENCES provided.

        ELM_MATCH_VEC = ELMMatchSeqs()


        ELM_MATCH_VEC   A LIST corresponding to each BioSeq object provided.  Each element in the list is a
                        DICT where each element indicates a matched ELM location.

        """

        if len(self.ELMDict) == 0:
            print 'ELMParser must be run first!'
            return
        if len(self.Sequences) == 0:
            print 'LoadSeqs must be run first!'
            return


        if len(self.PreCalcELMMatchSeqs) == 0:
            #compile all of the re's
            if self.ELMreDict != None:
                ELMreDict={}
                for thisELM in self.ELMOrder:
                    ELMreDict[thisELM]=re.compile(self.ELMDict[thisELM],re.I)
                self.ELMreDict=ELMreDict
            
            SeqELMDict=[]
            for thisSeq in self.Sequences:
                thisELMSeq={}
                for thisELM in self.ELMOrder:
                    thisIndex=[];
                    #each search only finds one match, so repeat the search until no more are found
                    spot = self.ELMreDict[thisELM].search(thisSeq)
                    while spot != None:
                        thisIndex.append(spot.start())
                        spot = self.ELMreDict[thisELM].search(thisSeq,spot.start()+1)
                    thisELMSeq[thisELM]=thisIndex
                SeqELMDict.append(thisELMSeq)


            self.PreCalcELMMatchSeqs=SeqELMDict

        
        return self.PreCalcELMMatchSeqs
        
    def ELMHistGen(self):
        """
        ELMHistGen
            Generates a count-matrix of ELM matches at specific locations in the sequences.

        ELM_COUNT_MAT = ELMHistGen()

        ELM_COUNT_MAT   An [maxSeqLength x len(ELM_DICT)] each Row represent an ELM and each Column
                        is a seqLocation.

        """
        if len(self.ELMDict) == 0:
            print 'ELMParser must be run first!'
            return
        if len(self.Sequences) == 0:
            print 'LoadSeqs must be run first!'
            return

        if len(self.PreCalcELMMatchSeqs) == 0:
            self.ELMMatchSeqs()

        if len(self.PreCalchistArray) == 0:

            histArray=zeros((len(self.ELMOrder),self.MaxSize))

            counter=0
            for thisELM in self.ELMOrder:
                for i in xrange(len(self.PreCalcELMMatchSeqs)):
                    histArray[counter,self.PreCalcELMMatchSeqs[i][thisELM]] += 1
                counter+=1

            self.PreCalchistArray=histArray
        
        return self.PreCalchistArray

    def SimpleELMCall(self,SEQ_IND,WANTED_ELM):
        """
        SimpleELMCall
            Returns a simple Presence/Absence Call for the ELM.

        PACall = SimpleELMCall(SEQ_IND,WANTED_ELM)

        SEQ_IND     The index of the desired sequence
        WANTED_ELM  The desired ELM

        PACall      Returns 1 if the ELM is present in the sequence and 0 otherwise.

        """

        if SEQ_IND > len(self.Sequences) | SEQ_IND < 0:
            print 'Arguement SEQ_IND must be between 0 and len(self.Sequences)'
            raise IndexError
        
        if len(self.PreCalcELMMatchSeqs) == 0:
            self.ELMMatchSeqs()

        if len(self.PreCalcELMMatchSeqs[SEQ_IND][WANTED_ELM]) >= 1:
            return 1
        else:
            return 0



    def MultiELMCall(self,SEQ_IND,WANTED_ELM):
        """
        MultiELMCall
            Returns the binned Presence/Absence Call for the ELM.

        PACall = MultiELMCall(SEQ_IND,WANTED_ELM)

        SEQ_IND     The index of the desired sequence
        WANTED_ELM  The desired ELM

        PACall      Returns 1 if the ELM is present in the sequence and 0 otherwise.

        """
        if WANTED_ELM not in self.PreCalcELMBinDict:
            raise KeyError
            
        thisCall = []
        allBins = self.PreCalcELMBinDict[WANTED_ELM]
        for i in xrange(len(allBins)-1):
            if sum(self.ELMPosGen(SEQ_IND,allBins[i],allBins[i+1],WANTED_ELM)) > 0:
                thisCall.append(1)
            else:
                thisCall.append(0)
        return thisCall


        


    def ELMPosGen(self,SEQ_IND,POS_START,POS_STOP,WANTED_ELM):
        """
        ELMPosGen
            Generates a [len(ELM_MATCH_VEC) MAX_SIZE] matrix in which each position is a 1 if the ELM is present at that location
        and 0 otherwise.

            ELM_POS_DICT = ELMPosGen(POS_START,POS_STOP,WANTED_ELM)

        ELM_POS_DICT    A DICT which is 'keyed' by ELM name and each value is a matrix.


            ELM_POS_ARRAY = ELMPosGen(..., WANTED_ELM)
        Will return only the desired array.

        """

        #check to see if calculation can be performed
        if len(self.PreCalcELMMatchSeqs) == 0:
            self.ELMMatchSeqs()

        if POS_START < 0 | POS_STOP <= POS_START | POS_START > self.MaxSize | POS_STOP > self.MaxSize:
            print 'Arguements to POS_START >0 & POS_STOP > POS_START'
            raise IndexError

        if SEQ_IND > len(self.Sequences) | SEQ_IND < 0:
            print 'Arguement SEQ_IND must be between 0 and len(self.Sequences)'
            raise IndexError

        if not(WANTED_ELM in self.ELMOrder):
            print 'WANTED_ELM not found in ELMOrder'
            raise KeyError

        boolArray = zeros((1,POS_STOP-POS_START),'bool')
        for i in xrange(POS_START,POS_STOP):
            boolArray[0,i-POS_START]=i in self.PreCalcELMMatchSeqs[SEQ_IND][WANTED_ELM]
            
        return boolArray

    def GetELMIterator(self):
        """
        GetELMIterator
            Returns an iterator which will iterate over the items in the ELMDict.
        """
        return self.ELMIterator(self)
        

    class ELMIterator():
        def __init__(self,instance):
            self.__list=instance.ELMOrder
            self.__thisSpot=-1

            
        def __iter__(self):
            return self

        def next(self):
            self.__thisSpot += 1
            if self.__thisSpot < len(self.__list):
                return self.__list[self.__thisSpot]
            else:
                raise StopIteration
    
    def GetSeqIterator(self,WANTED_ELM=None,WANTED_RANGE=None):
        """
        GetSeqIterator
            Will allow for easy iteration and retrival of data in a Seq-Order manner.

        GetSeqIterator('ELMPosDict',WANTED_ELM,WANTED_RANGE)
            WANTED_RANGE        A tuple indicating the range of positions desired.  If left empty then the
                                entire sequence is provided.
            
            This creates an iterator which will return the ELMPos Array data in Sequence order. If WANTED_ELM
            is left empty then Data is returned as an Array in which each row represents one ELM.
        """
        if (WANTED_ELM == None) | (WANTED_ELM in self.ELMOrder):
            if WANTED_RANGE == None:
                WANTED_RANGE = (0,self.MaxSize-1)
            elif len(WANTED_RANGE) != 2:
                print 'WANTED_RANGE must be len == 2 or None'
                raise IndexError
            elif WANTED_RANGE[0] < 0 | WANTED_RANGE[1] <= WANTED_RANGE[0]:
                print 'WANTED_RANGE[0] > 0 and WANTED_RANGE[1] > WANTED_RANGE[0]'
                raise IndexError
            elif WANTED_RANGE[0] > self.MaxSize | WANTED_RANGE[1] > self.MaxSize:
                print 'WANTED_RANGE[0] < MaxSize and WANTED_RANGE[1] < MaxSize'
                raise IndexError
               
            return self.SeqIterator(self,WANTED_ELM,WANTED_RANGE)
        else:
            print 'Provided an unknown ELM key'
            raise KeyError
            

    class SeqIterator:
        def __init__(self,instance,WANTED_ELM,WANTED_RANGE):
            self.__WantedELM=WANTED_ELM
            self.__WantedRange=WANTED_RANGE
            self.__listSize=len(instance.Sequences)
            self.__thisSpot=-1
            self.__MaxSize=instance.MaxSize
            self.__ELMOrder=instance.ELMOrder
            self.ELMOrder=instance.ELMOrder
            self.PreCalcELMMatchSeqs=instance.PreCalcELMMatchSeqs
            self.ELMPosGen=instance.ELMPosGen

        def __iter__(self):
            return self

        def next(self):
            self.__thisSpot += 1
            if self.__thisSpot < self.__listSize:
                if self.__WantedELM != None:
                    return self.ELMPosGen(self.__thisSpot,self.__WantedRange[0],self.__WantedRange[1],self.__WantedELM)
                else:
                    OutArray=zeros((len(self.__ELMOrder),self.__WantedRange[1]-self.__WantedRange[0]))
                    for i in xrange(len(self.__ELMOrder)):
                        OutArray[i,:]=self.ELMPosGen(self.__thisSpot,self.__WantedRange[0],self.__WantedRange[1],self.__ELMOrder[i])
                    return OutArray
            else:
                raise StopIteration

    def __QueueWorker(self,WANTED_ELM):
        """
        A helper function to manage the Priority Queue and Threads
        """
        while True:
            try:
                thisItem = self.__WorkingQueue.get(True, 5)
                self.__WorkingQueue.task_done()
            except:
                print 'queue empty'
                break
            
            if not(self.__FoundCorrect):
                self.__MakeBreadth(WANTED_ELM,thisItem[0],thisItem[1],thisItem[2])
                

    def __MakeBreadth(self,WANTED_ELM,allBins,BreadthCounter,DepthCounter):
        """
        __MakeBreadth

        Performs a breadth-first search.  All nodes are added to the Queue and processes in order.


        """
        print (BreadthCounter,DepthCounter,self.__EvaluateAllBins(WANTED_ELM,allBins),self.__WorkingQueue.qsize(),WANTED_ELM)
        if DepthCounter > self.MaxDepth:
            return

        if (BreadthCounter > self.MaxBreadth) & (allBins in self.__DepthChecked):
            return
            
        LocalQueue = []
        LocalVals = []
            
        for k in xrange(len(allBins)-1):
            thisBinSet = deepcopy(allBins)
            newBinEdge=thisBinSet[k]+int((thisBinSet[k+1]-thisBinSet[k])/2)
            thisBinSet.insert(k+1,newBinEdge)
            LocalQueue.append(thisBinSet)
            LocalVals.append(self.__EvaluateAllBins(WANTED_ELM,thisBinSet))

        for k in xrange(1,len(allBins)-1):
            thisBinSet = deepcopy(allBins)
            thisBinSet.pop(k)
            LocalQueue.append(thisBinSet)
            LocalVals.append(self.__EvaluateAllBins(WANTED_ELM,thisBinSet))

        for k in xrange(1,len(allBins)-1):
            thisBinSet = deepcopy(allBins)
            newBinMove = int((thisBinSet[k+1]-thisBinSet[k])/2)
            LocalQueue.append(thisBinSet)
            LocalVals.append(self.__EvaluateAllBins(WANTED_ELM,thisBinSet))

        for k in xrange(1,len(allBins)-1):
            thisBinSet = deepcopy(allBins)
            newBinMove=int((thisBinSet[k]-thisBinSet[k-1])/2)
            thisBinSet[k]=thisBinSet[k]-newBinMove
            LocalQueue.append(thisBinSet)
            LocalVals.append(self.__EvaluateAllBins(WANTED_ELM,thisBinSet))

        BestInds=argsort(LocalVals)


        if BreadthCounter < self.MaxBreadth:
            for k in xrange(min(self.MaxWidth,len(BestInds))-1,0,-1):
                print 'inside here'
                self.__WorkingQueue.put(((LocalQueue[BestInds[k]],BreadthCounter+1,DepthCounter),LocalVals[BestInds[k]]))
        else:
            if (LocalVals[BestInds[0]] < self.__EvaluateAllBins(WANTED_ELM,allBins)):
                self.__DepthChecked.append(allBins)
                self.__WorkingQueue.put(((LocalQueue.pop(BestInds[0]),BreadthCounter+1,DepthCounter+1),LocalVals[BestInds[0]]))


    

    def __EvaluateBin(self,WANTED_ELM,binEdges):
        """
        Evaluates the values of the single bin provided.
        """
        thisBinVal=[0,0,0]
        for thisSeq in self.GetSeqIterator(WANTED_ELM,(binEdges[0],binEdges[1])):
            CurrentBin=sum(thisSeq)
            if CurrentBin == 0:
                thisBinVal[0]+=1
            elif CurrentBin == 1:
                thisBinVal[1]+=1
            else:
                thisBinVal[2]+=1
        return (thisBinVal[0], thisBinVal[1], thisBinVal[2] )


    def __EvaluateAllBins(self,WANTED_ELM,allBins):
        """
        Evaluates the values of the bins provided.  Uses a dictionary to avoid re-calculating values wherever possible
        """
        if not(tuple(allBins) in self.__allBinDict ):
            TotalEmpty=0
            TotalCorrect=0
            TotalWrong=0
            for i in xrange(len(allBins)-1):
                if (allBins[i],allBins[i+1]) not in self.__singlebinDict:
                    self.__singlebinDict[(allBins[i],allBins[i+1])] = self.__EvaluateBin(WANTED_ELM,(allBins[i],allBins[i+1]))
                
                TotalEmpty += self.__singlebinDict[(allBins[i],allBins[i+1])][0]
                TotalCorrect += self.__singlebinDict[(allBins[i],allBins[i+1])][1]
                TotalWrong += self.__singlebinDict[(allBins[i],allBins[i+1])][2]
            funVal = (.25*TotalEmpty + TotalWrong) / (TotalCorrect + 1)
            self.__allBinDict[tuple(allBins)]=(TotalEmpty,TotalCorrect,TotalWrong,funVal)
            if TotalWrong == 0:
                self.__FoundCorrect = True
        return self.__allBinDict[tuple(allBins)][3]

    def CalculateBins(self,WANTED_ELM):
        if len(self.PreCalcELMMatchSeqs) == 0:
            self.ELMMatchSeqs()

        self.__singlebinDict={}
        self.__allBinDict={}
        self.__DepthChecked=[]
        self.__FoundCorrect = False
        
        self.__WorkingQueue = PriorityQueue(-1)
        theseBins = [0, int(self.MaxSize/2), self.MaxSize]

        self.__WorkingQueue.put(((theseBins,0,0),0))

        self.__QueueWorker(WANTED_ELM)

        currentMin = 50000
        currentMinInd = 0
        for i in self.__allBinDict:
            if self.__allBinDict[i][3] < currentMin:
                currentMinInd = i
                currentMin = self.__allBinDict[currentMinInd][3]

        self.PreCalcELMBinDict[WANTED_ELM] = currentMinInd
