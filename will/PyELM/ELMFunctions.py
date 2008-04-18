import os

def Parser(DIRECTORY="C:\\Documents and Settings\\Will\\My Documents\\PyELM\\ELM_RAW_DOWNLOAD\\"):
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
