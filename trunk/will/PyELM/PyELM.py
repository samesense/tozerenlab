import os
from numpy import *
import re
from pylab import *

print 'Actually running'


def Parser(DIRECTORY="C:\Documents and Settings\Will\My Documents\ELM_Motif_finder\ELM_RAW_DOWNLOAD\\"):
    """
    Parser
        Parses a directory of HTML files to extract the ELM names and Regular Expressions and returns a DICT.

        ELM_LIST = Parser(DIRECTORY)

        DIRECTORY       The directory containing the HTML files downloaded from the ELM database.
                        If left empty then C:\Documents and Settings\Will\My Documents\PyELM\ELM_RAW_DOWNLOAD\

        ELM_LIST        A DICT containing all of the parsed ELMs.
    """


    fileList=os.listdir(DIRECTORY)
    allRegExps={}
    
    for i in xrange(1,len(fileList)-1):
        
        thisFile=open(DIRECTORY+fileList[i],'r')
        ELMname=fileList[i].rpartition('.')
        currentLine=''
        
        while currentLine.find('Pattern:') == -1:
            currentLine=thisFile.readline()
        
        currentLine=thisFile.readline()
        
        allRegExps[ELMname[0]]=currentLine[currentLine.find('>')+1:currentLine.find('/')-2]
        thisFile.close()
    
    return allRegExps


    fileList=os.listdir(DIRECTORY)
    allRegExps={}

    for i in xrange(1,len(fileList)-1):
        
        thisFile=open(DIRECTORY+fileList[i],'r')
        ELMname=fileList[i].rpartition('.')
        currentLine=''
        
        while currentLine.find('Pattern:') == -1:
            currentLine=thisFile.readline()
        
        currentLine=thisFile.readline()
        
        allRegExps[ELMname[0]]=currentLine[currentLine.find('>')+1:currentLine.find('/')-3]
        thisFile.close()
    
    return allRegExps

def ELMMatchSeqs(SEQUENCES,ELM_DICT):
    """
    ELMMatchSeqs
        Searches for the regular expressions of ELM_DICT within the SEQUENCES provided.

    ELM_MATCH_VEC = ELMMatchSeqs(SEQUENCES,ELM_DICT)

    SEQUENCES       A LIST of AA sequences.

    ELM_DICT        A DICT of ELM motifs as created by Parser()

    ELM_MATCH_VEC   A LIST corresponding to each BioSeq object provided.  Each element in the list is a
                    DICT where each element indicates a matched ELM location.

    """

    #compile all of the re's
    ELMreDict={}
    for thisELM in ELM_DICT:
        ELMreDict[thisELM]=re.compile(ELM_DICT[thisELM],re.I)

    
    SeqELMDict=[]
    for thisSeq in SEQUENCES:
        thisELMSeq={}
        for thisELM in ELMreDict:
            thisIndex=[];
            #each search only finds one match, so repeat the search until no more are found
            spot = ELMreDict[thisELM].search(thisSeq)
            while spot != None:
                thisIndex.append(spot.start())
                spot = ELMreDict[thisELM].search(thisSeq,spot.start()+1)
            thisELMSeq[thisELM]=thisIndex
        SeqELMDict.append(thisELMSeq)

    return SeqELMDict
    
def ELMHistGen(ELM_MATCH_VEC,ELM_DICT,MAX_SIZE):
    """
    ELMHistGen
        Generates a count-matrix of ELM matches at specific locations in the sequences.

    ELM_COUNT_MAT = ELMHistGen(ELM_MATCH_VEC,ELM_DICT,MAX_SIZE)

    ELM_MATCH_VEC   An output produced by ELMMatchSeqs.  This is an array where each element represnts a
                    sequence and each value is a DICT of ELMs and thier matched locations.

    ELM_DICT        An output of Parser().  This is a DICT of each ELM and its regexp.

    ELM_COUNT_MAT   An [maxSeqLength x len(ELM_DICT)] each Row represent an ELM and each Column
                    is a seqLocation.

    """

    histArray=zeros((len(ELM_DICT),MAX_SIZE))

    counter=0
    for thisELM in ELM_DICT:
        for i in xrange(len(ELM_MATCH_VEC)):
            histArray[counter,ELM_MATCH_VEC[i][thisELM]] += 1
        counter+=1
    
    return histArray

def ELMPosGen(ELM_MATCH_VEC,ELM_DICT,MAX_SIZE,WANTED_ELM=None):
    """
    ELMPosGen
        Generates a [len(ELM_MATCH_VEC) MAX_SIZE] matrix in which each position is a 1 if the ELM is present at that location
    and 0 otherwise.

        ELM_POS_DICT = ELMPosGen(ELM_MATCH_VEC,ELM_DICT,MAX_SIZE)

    ELM_MATCH_VEC   An output produced by ELMMatchSeqs.  This is an array where each element represnts a
                    sequence and each value is a DICT of ELMs and thier matched locations.

    ELM_DICT        An output of Parser().  This is a DICT of each ELM and its regexp.

    ELM_POS_DICT    A DICT which is 'keyed' by ELM name and each value is a matrix.


        ELM_POS_ARRAY = ELMPosGen(..., WANTED_ELM)
    Will return only the desired array.

    """
    

    if WANTED_ELM != None:
        boolArray = zeros((len(ELM_MATCH_VEC),MAX_SIZE))
        for i in xrange(len(ELM_MATCH_VEC)):
            boolArray[i,ELM_MATCH_VEC[i][WANTED_ELM]]=1
        return boolArray

    ELM_POS_DICT={}
    for thisELM in ELM_DICT:
        boolArray = zeros((len(ELM_MATCH_VEC),MAX_SIZE))
        for i in xrange(len(ELM_MATCH_VEC)):
            boolArray[i,ELM_MATCH_VEC[i][thisELM]]=1
        ELM_POS_DICT[thisELM] = boolArray

    return ELM_POS_DICT

    

    



















    
