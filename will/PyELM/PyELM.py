import os
import Bio
import numpy
import re

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

    SEQUENCES       A LIST or ITER of BioSeq objects of Amino Acid sequences

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
            spot = ELMreDict[thisELM].search(thisSeq.seq.tostring())
            while spot != None:
                thisIndex.append(spot.start())
                spot = ELMreDict[thisELM].search(thisSeq.seq.tostring(),spot.start()+1)
            thisELMSeq[thisELM]=thisIndex
        SeqELMDict.append(thisELMSeq)

    return SeqELMDict
    
