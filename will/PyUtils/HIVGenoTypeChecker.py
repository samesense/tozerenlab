import re
import os
from PyMozilla import *
from Bio import SeqIO
import time


def GetHIVGenotype(INPUT_SEQ):
    """
    GetHIVGenotype(INPUT_SEQ)
        INPUT_SEQ       A string of nucleotides from an HIV sequence

    Returns:
        Tuple (GI_NAMES,NUMERIC_DATA,NUM_WINDOWS)

        GI_NAMES        An ordered (by most likely subtype) list of HIV-1 subtypes
        NUMERIC_DATA    The windowed BLAST scores
        NUM_WINDOWS     The number of windows that the sequence was divided into

        If there was an error then it returns an empty tuple ()
    """

    try:    
        MozEmu = MozillaEmulator(cacher=None,trycount=10)
        thisPost = 'BLAST_DATABASE=1&QUERY_SEQUENCE='+INPUT_SEQ
        retData = MozEmu.download('http://www.ncbi.nlm.nih.gov/projects/genotyping/genotype.cgi',postdata=thisPost)

    except:
        return ()
   
    try:
        GI_NAMES=re.search('name="GI_POINTER" value="(.*)"',retData).group(1)
        NUMERIC_DATA=re.search('name="SCORE_ARRAY" value="(.*)"',retData).group(1)
        NUM_WINDOWS=re.search('name="NUMBER_OF_WINDOWS" value="(.*)"',retData).group(1)
        outputData = (GI_NAMES,NUMERIC_DATA,NUM_WINDOWS)
        print 'Success'
    except:
        print 'Failure'
        outputData = ()
        
    return outputData

def WriteData(SIMPLE_FHANDLE,FULL_FHANDLE,NAME,HIV_GENO_DATA):
    """
    WriteData(SIMPLE_FHANDLE,FULL_FHANDLE,HIV_GENO_DATA)
            SIMPLE_FHANDLE      A file handle to the Simple Output Data file
            FULL_FHANDLE        A file handle to the Full Output Data file
            NAME                A header name
            HIV_GENO_DATA       The tuple output of GetHIVGenotype

    Returns:
        None

    Simple Output
        HEADER,Likely Genotype \n

    Full Output
        Header \n
        GI_NAMES \n
        NUMERIC_DATA \n
        NUM_WINDOWS \n

    """

    BestSub = re.split(',',HIV_GENO_DATA[0])[0]

    SIMPLE_FHANDLE.write(NAME+','+BestSub+'\n')

    FULL_FHANDLE.write(NAME+'\n')
    FULL_FHANDLE.write(HIV_GENO_DATA[0]+'\n')
    FULL_FHANDLE.write(HIV_GENO_DATA[1]+'\n')
    FULL_FHANDLE.write(HIV_GENO_DATA[2]+'\n')

    


def CheckAllGenoTypes(FASTA_DIR,FULL_FILE,SIMPLE_FILE,MISSED_FILE):
    FullFile = open(FULL_FILE,'w')
    SimpleFile = open(SIMPLE_FILE,'w')
    MissedFile = open(MISSED_FILE,'w')
    counter = 0

    for thisSeq in FastaDirIter(FASTA_DIR):
        counter+=1
        
        print counter
        GenoData = GetHIVGenotype(thisSeq.seq.tostring())
        
        if len(GenoData) == 3:
            WriteData(SimpleFile,FullFile,thisSeq.description,GenoData)
        else:
            MissedFile.write(thisSeq.description+'\n'+thisSeq.seq.tostring())
        #time.sleep(30)


    FullFile.close()
    SimpleFile.close()
    MissedFile.close()
    
            
	



    



class FastaDirIter:
    def __init__(self,DIRECTORY):
        self.Files=os.listdir(DIRECTORY)
        if DIRECTORY[len(DIRECTORY)-1] != '\\' :
            self.Directory=DIRECTORY+'\\'
        else:
            self.Directory=DIRECTORY
            
        self.CurrentHandle = open(self.Directory+self.Files.pop(),'r')
        self.CurrentIter = SeqIO.parse(self.CurrentHandle,'fasta')
        

    def __iter__(self):
        return self

    def next(self):
        try:
            return self.CurrentIter.next()
        except StopIteration:
            self.CurrentHandle.close()
            if len(self.Files)>0:
                self.CurrentHandle = open(self.Directory+self.Files.pop(),'r')
                self.CurrentIter = SeqIO.parse(self.CurrentHandle,'fasta')
                return self.CurrentIter.next()
            else:
                raise StopIteration
        
        

    



    
