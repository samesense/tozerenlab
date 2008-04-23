import os
from numpy import *
import re
from pylab import *

print 'Actually running'

class PyELM:

    def __init__(self):

        self.ELMDict={}
        self.ELMOrder=[]
        self.Sequences=[]
        self.ELMreDict={}
        self.MaxSize=0    


        self.PreCalchistArray=[]        
        self.PreCalcELMMatchSeqs=[]

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


    

    def ELMMatchSeqs(self,SEQUENCES=None, ELM_DICT=None):
        """
        ELMMatchSeqs
            Searches for the regular expressions of ELM_DICT within the SEQUENCES provided.

        ELM_MATCH_VEC = ELMMatchSeqs(SEQUENCES,ELM_DICT)


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
            elif WANTED_RANGE[0] < 0 | WANTED_RANGE[1] < WANTED_RANGE[0]:
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
