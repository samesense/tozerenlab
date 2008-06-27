"""
HIVGenoTypeChecker
    A module that will interact with the NCBI HIV Genotyping tool at
    'http://www.ncbi.nlm.nih.gov/projects/genotyping/genotype.cgi' and determine
    the contributing HIV-1 subtypes with a set of sequences.
"""

import re
import os
import PyMozilla 
from Bio import SeqIO
from urllib2 import URLError



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

    ncbi_cgi = 'http://www.ncbi.nlm.nih.gov/projects/genotyping/genotype.cgi'
    try:    
        mozEmu = PyMozilla.MozillaEmulator(cacher = None, trycount = 10)
        thisPost = 'BLAST_DATABASE=1&QUERY_SEQUENCE=' + INPUT_SEQ
        
        retData = mozEmu.download(ncbi_cgi, postdata = thisPost)

    except URLError:
        return ()
   
    try:
        gi_names = re.search('name="GI_POINTER" value="(.*)"', retData).group(1)

        numeric_data = re.search('name="SCORE_ARRAY" value="(.*)"',
                               retData).group(1)

        num_windows = re.search('name="NUMBER_OF_WINDOWS" value="(.*)"',
                              retData).group(1)

        outputData = (gi_names, numeric_data, num_windows)
        print 'Success'
    except AttributeError:
        print 'Failure'
        outputData = ()
        
    return outputData

def WriteData(SIMPLE_FHANDLE, FULL_FHANDLE, NAME, HIV_GENO_DATA):
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

    bestSub = re.split(',', HIV_GENO_DATA[0])[0]

    SIMPLE_FHANDLE.write(NAME + ',' + bestSub + '\n')

    FULL_FHANDLE.write(NAME + '\n')
    FULL_FHANDLE.write(HIV_GENO_DATA[0] + '\n')
    FULL_FHANDLE.write(HIV_GENO_DATA[1] + '\n')
    FULL_FHANDLE.write(HIV_GENO_DATA[2] + '\n')

    


def CheckAllGenoTypes(FASTA_DIR, FULL_FILE, SIMPLE_FILE, MISSED_FILE):
    """
    CheckAllGenoTypes
        Does all of the required actions for checking HIV-1 Genotypes.

    CheckAllGenoType(FASTA_DIR, FULL_FILE, SIMPLE_FILE, MISSED_FILE)


    """
    fullFileHandle = open(FULL_FILE, 'w')
    simpleFileHandle = open(SIMPLE_FILE, 'w')
    missedFileHandle = open(MISSED_FILE, 'w')
    counter = 0

    for thisSeq in FastaDirIter(FASTA_DIR):
        counter += 1
        
        print counter
        genoData = GetHIVGenotype(thisSeq.seq.tostring())
        
        if len(genoData) == 3:
            WriteData(simpleFileHandle, fullFileHandle,
                      thisSeq.description, genoData)
        else:
            missedFileHandle.write(thisSeq.description +
                                   '\n' + thisSeq.seq.tostring())
        #time.sleep(30)


    fullFileHandle.close()
    simpleFileHandle.close()
    missedFileHandle.close()
    
            
	



    



class FastaDirIter:
    """
    FastaDirIter
        An iterator class which steps through multiple fasta files in one
    directory and then passes the SeqRecord back at each iteration.

    """
    def __init__(self, DIRECTORY):
        """
        __init__(self,DIRECTORY)

        DIRECTORY defines the path to a directory of fasta files.
        """
        self.files = os.listdir(DIRECTORY)
        if DIRECTORY[len(DIRECTORY)-1] != os.sep :
            self.directory = DIRECTORY + os.sep
        else:
            self.directory = DIRECTORY
            
        self.current_handle = open(self.directory + self.files.pop(), 'r')
        self.current_iter = SeqIO.parse(self.current_handle, 'fasta')
        

    def __iter__(self):
        """
        __iter__(self)
            Returns the class for iterations 

        """
        return self

    def next(self):
        """
        next(self)
            Returns the next record in the sequence.

        """
        try:
            return self.current_iter.next()
        except StopIteration:
            self.current_handle.close()
            if len(self.files)>0:
                self.current_handle = open(self.directory + self.files.pop(),
                                          'r')
                
                self.current_iter = SeqIO.parse(self.current_handle, 'fasta')
                return self.current_iter.next()
            else:
                raise StopIteration
        
        

    



    
