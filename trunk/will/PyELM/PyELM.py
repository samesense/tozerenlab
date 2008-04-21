import os
import Bio
import numpy
import re

print 'Actually running'


def Parser(DIRECTORY="C:\\Documents and Settings\\Will\\My Documents\\PyELM\\ELM_RAW_DOWNLOAD\\"):
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
        
        allRegExps[ELMname[0]]=currentLine[currentLine.find('>')+1:currentLine.find('/')-3]
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
