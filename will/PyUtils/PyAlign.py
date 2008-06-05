# -*- coding: cp1252 -*-
#from queue import *
import tempfile
import os
import subprocess
import re

class ClustalInterface():
    def __init__(self,clustalPath='C:\\clustalw\\',InputSeqs=[],ParamPairs=[]):
        self.clustalPath = clustalPath
        self.InputSeqs = InputSeqs
        self.FastaName = None
        self.OutputName = ''
        self.ParamPairs=[]
        self.alignedSeqs=[];

        self.FastaWritten = False
        self.ProcessRun = False
        self.CleanedUp = False
        

        self.ParamDict = {}
        self.MakeParamDict()


    def WriteFASTA(self,InputSeqs=None):
        """
        WriteFASTA(InputSeqs=None)
            Uses NamedTemporaryFile to generate a unique filename and writes the sequences in FASTA format.  If InputSeqs is not empty then they are added to the class.


        """

        
        newNamedFile = tempfile.NamedTemporaryFile(prefix='temp_',suffix='.fasta',dir=self.clustalPath)

        self.FastaName = newNamedFile.name
        newNamedFile.close()
        
        FastaHandle = open(self.FastaName,'wt')

        for i in xrange(len(self.InputSeqs)):
            FastaHandle.write('>Seq' + str(i) + '\n')
            FastaHandle.write(self.InputSeqs[i] + '\n')
            self.alignedSeqs.append('')

        FastaHandle.close()

        self.FastaWritten = True
        

    def MakeParamDict(self):
        """
        MakeParamDict()
            Creates the parameter dictionary for the ClustalW interface.  Changing these can be modified to change the functionality of the ClustalW interface.

        ***General settings:****
        'INTERACTIVE':      read command line, then enter normal interactive menus
        'QUICKTREE':        use FAST algorithm for the alignment guide tree
        'NEGATIVE':         protein alignment with negative values in matrix
        'OUTFILE':          sequence alignment file name
    
        'OUTORDER':         [INPUT] or ALIGNED
        'CASE':             [LOWER] or UPPER (for GDE output only)
        'SEQNOS':           [OFF] or ON (for Clustal output only)

        ***Fast Pairwise Alignments:***
        'KTUPLE':           word size [1]
        'TOPDIAGS':         number of best diags. [5]
        'WINDOW':           window around best diags. [5]
        'PAIRGAP':          gap penalty [3]
        'SCORE':            PERCENT or [ABSOLUTE]

        ***Slow Pairwise Alignments:***
        'PWMATRIX':         Protein weight matrix=BLOSUM, PAM, [GONNET], ID or filename
        'PWDNAMATRIX':      DNA weight matrix=[IUB], CLUSTALW or filename²
        'PWGAPOPEN':        gap opening penalty [10]
        'PWGAPEXT':         gap extension penalty [0.1]

        ***Multiple Alignments:***
        'NEWTREE':          file for new guide tree
        'MATRIX':           Protein weight matrix=BLOSUM, PAM, [GONNET], ID or filename
        'DNAMATRIX':        DNA weight matrix=[IUB], CLUSTALW or filename
        'GAPOPEN':          gap opening penalty [10]
        'GAPEXT':           gap extension penalty [0.2]
        'ENDGAPS':          no end gap separation pen.
        'GAPDIST':          gap separation pen. range [4]
        'NOPGAP':           residue-specific gaps off
        'NOHGAP':           hydrophilic gaps off
        'HGAPRESIDUES':     list hydrophilic res. ['GPSNDQEKR']
        'MAXDIV':           % ident. for delay [30]
        'TYPE':             PROTEIN or DNA [AUTO]
            

        """

        self.ParamDict['INFILE']=None               #input sequences

        
        #***General settings:****
        self.ParamDict['INTERACTIVE']=None          #:read command line, then enter normal interactive menus
        self.ParamDict['QUICKTREE']=None            #:use FAST algorithm for the alignment guide tree
        self.ParamDict['NEGATIVE']=None             #:protein alignment with negative values in matrix
        self.ParamDict['OUTFILE']=None              #:sequence alignment file name
    
        self.ParamDict['OUTORDER']=None             #:[INPUT] or ALIGNED
        self.ParamDict['CASE']=None                 #:[LOWER] or UPPER (for GDE output only)
        self.ParamDict['SEQNOS']=None               #:[OFF] or ON (for Clustal output only)

        #***Fast Pairwise Alignments:***
        self.ParamDict['KTUPLE']=None               #:word size [1]
        self.ParamDict['TOPDIAGS']=None             #:number of best diags. [5]
        self.ParamDict['WINDOW']=None               #:window around best diags. [5]
        self.ParamDict['PAIRGAP']=None              #:gap penalty [3]
        self.ParamDict['SCORE']=None                #:PERCENT or [ABSOLUTE]

        #***Slow Pairwise Alignments:***
        self.ParamDict['PWMATRIX']=None             #:Protein weight matrix=BLOSUM, PAM, [GONNET], ID or filename
        self.ParamDict['PWDNAMATRIX']=None          #:DNA weight matrix=[IUB], CLUSTALW or filename²
        self.ParamDict['PWGAPOPEN']=None            #:gap opening penalty [10]
        self.ParamDict['PWGAPEXT']=None             #:gap extension penalty [0.1]

        #***Multiple Alignments:***
        self.ParamDict['NEWTREE']=None              #:file for new guide tree
        self.ParamDict['MATRIX']=None               #:Protein weight matrix=BLOSUM, PAM, [GONNET], ID or filename
        self.ParamDict['DNAMATRIX']=None            #:DNA weight matrix=[IUB], CLUSTALW or filename
        self.ParamDict['GAPOPEN']=None              #:gap opening penalty [10]
        self.ParamDict['GAPEXT']=None               #:gap extension penalty [0.2]
        self.ParamDict['ENDGAPS']=None              #:no end gap separation pen.
        self.ParamDict['GAPDIST']=None              #:gap separation pen. range [4]
        self.ParamDict['NOPGAP']=None               #:residue-specific gaps off
        self.ParamDict['NOHGAP']=None               #:hydrophilic gaps off
        self.ParamDict['HGAPRESIDUES']=None         #:list hydrophilic res. ['GPSNDQEKR']
        self.ParamDict['MAXDIV']=None               #:% ident. for delay [30]
        self.ParamDict['TYPE']=None                 #:PROTEIN or DNA [AUTO]


    def MakeCommand(self):
        """
        MakeCommand()
            Creates the command required to run the ClustalW.  It will append any parameters from self.ParamDict which are not None

        """
        command = 'clustalw ' + self.FastaName 

        for thisPair in self.ParamPairs:
            self.ParamDict[thisPair[0]]=str(thisPair[1])

        for thisParam in self.ParamDict:
            if self.ParamDict[thisParam] != None:
                command += ' -' + thisParam + '=' + self.ParamDict

        return command

    def DoAlignment(self):
        """
        DoAlignment()
            Uses the subprocess.Popen to run the ClustalW command and .wait() to block the interpeter until execution finishes.
        """

        processVar = subprocess.Popen(self.MakeCommand(),shell=True)
        processVar.wait()

        self.ProcessRun=True
        

    def ReadALN(self):
        """
        ReadALN()
            Reads the .aln file that resulted from the completion of DoAlignment().  It places the data into the self.alignedSeqs in the same order that they were provided to the class.
        """
        
        alnName = self.FastaName[0:-6] + '.aln'

        reFilter = re.compile('Seq(\d*) *(.*)')
        #returns the 'n' of Seqn and the subsequent alignment string

        alnHandle = open(alnName,'r+t')
        for line in alnHandle:
            reOutput = reFilter.findall(line)
            if len(reOutput) == 1:
                self.alignedSeqs[int(reOutput[0][0])]+=reOutput[0][1]


        alnHandle.close()

    def CleanUpFiles(self,fasta=True,aln=True,dnd=True):
        """
        CleanUpFiles(fasta=True,aln=True,dnd=True)
            Cleans up the .fasta; .aln; and .dnd files that are left over after the alignment finishes.

        """
        if fasta==True:
            try:
                os.remove(self.FastaName)
            except:
                print 'Fasta Already Removed'
        if aln:
            try:
                os.remove(self.FastaName[0:-6] + '.aln')
            except:
                print 'ALN Already Removed'
        if dnd:
            try:
                os.remove(self.FastaName[0:-6] + '.dnd')
            except:
                print 'DND Already Removed'
        self.CleanedUp=fasta&aln&dnd
        
    def AlignSeqs(self,InputSeqs=None,ParamPairs=None):
        """
        AlignSeqs(InputSeqs=None,ParamPairs=None)
            A gateway function which takes a set of Input Sequences and Parameter Pairs and performs the Writing, Alignment, and Reading in one Call.

        """

        if InputSeqs != None:
            self.CleanUpFiles(aln=False,dnd=False)
            self.FastaWritten = False
            self.InputSeqs=InputSeqs

        if ParamPairs != None:
            self.ParamPairs = ParamPairs

        if not(self.FastaWritten):
            self.WriteFASTA()
            
        self.DoAlignment()
        self.ReadALN()
        self.CleanUpFiles()
        
        


