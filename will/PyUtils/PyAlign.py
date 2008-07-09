# -*- coding: cp1252 -*-
"""
PyAlign
    Uses the subprocesses module to create an interface between python and the
    ClustalW command-line tool
"""

#from queue import *
import tempfile
import os
import subprocess
import re


def AlignSeqs(INPUT_SEQS, PARAMS = None):
    """
    AlignSeqs
        A shortcut to aligning using the default clustalw parameters.
    """
    cont = ClustalInterface(INPUT_SEQS)

    return cont.AlignSeqs()




class ClustalInterface():
    """
    ClustalInterface
        An interface for doing alignments with ClustalW.  The Class interface
        can write fasta-files, call clustalW and read the .aln file that results.

    """

    
    def __init__(self, CLUSTALPATH = 'C:\\clustalw\\', INPUT_SEQS = None,
                 PARAMPAIRS = None):
        """ 
        __init__(CLUSTALPATH = 'C:\\clustalw\\', INPUT_SEQS = None,
                 PARAMPAIRS = None)
        """
        
        self.clustalpath = CLUSTALPATH
        self.inputseqs = INPUT_SEQS
        self.fastaname = None
        self.outputname = ''
        self.parampairs = PARAMPAIRS
        self.alignedseqs = []

        self.fastawritten = False
        self.processrun = False
        self.cleanedup = False
        

        self.paramdict = {}
        self.MakeParamDict()


    def WriteFASTA(self, INPUT_SEQS = None):
        """
        WriteFASTA(inputSeqs = None)
            Uses NamedTemporaryFile to generate a unique filename and writes
        the sequences in FASTA format.  If InputSeqs is not empty then they
        are added to the class.


        """

        
        newNamedFile = tempfile.NamedTemporaryFile(prefix = 'temp_',
                                                   suffix = '.fasta',
                                                   dir = self.clustalpath)

        self.fastaname = newNamedFile.name
        newNamedFile.close()
        
        fastaHandle = open(self.fastaname, 'wt')

        for i in xrange(len(self.inputseqs)):
            fastaHandle.write('>Seq' + str(i) + '\n')
            fastaHandle.write(self.inputseqs[i] + '\n')
            self.alignedseqs.append('')

        fastaHandle.close()

        self.fastawritten = True
        

    def MakeParamDict(self):
        """
        MakeParamDict()
            Creates the parameter dictionary for the ClustalW interface.
        Changing these can be modified to change the functionality of the
        ClustalW interface.

        ***General settings:****
        'INTERACTIVE':      read command line, then enter normal
                            interactive menus
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
        'PWMATRIX':         Protein weight matrix=BLOSUM, PAM, [GONNET], ID or
                            filename
        'PWDNAMATRIX':      DNA weight matrix=[IUB], CLUSTALW or filename²
        'PWGAPOPEN':        gap opening penalty [10]
        'PWGAPEXT':         gap extension penalty [0.1]

        ***Multiple Alignments:***
        'NEWTREE':          file for new guide tree
        'MATRIX':           Protein weight matrix=BLOSUM, PAM, [GONNET], ID or
                            filename
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

        self.paramdict['INFILE'] = None               #input sequences

        
        #***General settings:****
        #:read command line, then enter normal interactive menus
        self.paramdict['INTERACTIVE'] = None
        #:use FAST algorithm for the alignment guide tree
        self.paramdict['QUICKTREE'] = None            
        #:protein alignment with negative values in matrix
        self.paramdict['NEGATIVE'] = None             
        #:sequence alignment file name
        self.paramdict['OUTFILE'] = None              
    
        #:[INPUT] or ALIGNED
        self.paramdict['OUTORDER'] = None
        #:[LOWER] or UPPER (for GDE output only)
        self.paramdict['CASE'] = None
        #:[OFF] or ON (for Clustal output only)
        self.paramdict['SEQNOS'] = None               

        #***Fast Pairwise Alignments:***
        #:word size [1]
        self.paramdict['KTUPLE'] = None
        #:number of best diags. [5]
        self.paramdict['TOPDIAGS'] = None
        #:window around best diags. [5]
        self.paramdict['WINDOW'] = None
        #:gap penalty [3]
        self.paramdict['PAIRGAP'] = None
        #:PERCENT or [ABSOLUTE]
        self.paramdict['SCORE'] = None                

        #***Slow Pairwise Alignments:***
        #:Protein weight matrix=BLOSUM, PAM, [GONNET], ID or filename
        self.paramdict['PWMATRIX'] = None
        #:DNA weight matrix=[IUB], CLUSTALW or filename²
        self.paramdict['PWDNAMATRIX'] = None
        #:gap opening penalty [10]
        self.paramdict['PWGAPOPEN'] = None
        #:gap extension penalty [0.1]
        self.paramdict['PWGAPEXT'] = None             

        #***Multiple Alignments:***
        #:file for new guide tree
        self.paramdict['NEWTREE'] = None
        #:Protein weight matrix=BLOSUM, PAM, [GONNET], ID or filename
        self.paramdict['MATRIX'] = None
        #:DNA weight matrix=[IUB], CLUSTALW or filename
        self.paramdict['DNAMATRIX'] = None
        #:gap opening penalty [10]
        self.paramdict['GAPOPEN'] = None
        #:gap extension penalty [0.2]
        self.paramdict['GAPEXT'] = None
        #:no end gap separation pen.
        self.paramdict['ENDGAPS'] = None
        #:gap separation pen. range [4]
        self.paramdict['GAPDIST'] = None
        #:residue-specific gaps off
        self.paramdict['NOPGAP'] = None
        #:hydrophilic gaps off
        self.paramdict['NOHGAP'] = None
        #:list hydrophilic res. ['GPSNDQEKR']
        self.paramdict['HGAPRESIDUES'] = None
        #:% ident. for delay [30]
        self.paramdict['MAXDIV'] = None
        #:PROTEIN or DNA [AUTO]
        self.paramdict['TYPE'] = None                 


    def MakeCommand(self):
        """
        MakeCommand()
            Creates the command required to run the ClustalW.  It will append
        any parameters from self.paramdict which are not None

        """
        command = 'clustalw ' + self.fastaname

        if self.parampairs != None:
            for thisPair in self.parampairs:
                self.paramdict[thisPair[0]] = str(thisPair[1])

        for thisParam in self.paramdict:
            if self.paramdict[thisParam] != None:
                command += ' -' + thisParam + '=' + self.paramdict

        return command

    def DoAlignment(self):
        """
        DoAlignment()
            Uses the subprocess.Popen to run the ClustalW command and .wait()
        to block the interpeter until execution finishes.
        """

        processVar = subprocess.Popen(self.MakeCommand(), shell = True)
        processVar.wait()

        self.processrun = True
        

    def ReadALN(self):
        """
        ReadALN()
            Reads the .aln file that resulted from the completion of
        DoAlignment().  It places the data into the self.alignedseqs in the
        same order that they were provided to the class.
        """
        
        alnName = self.fastaname[0:-6] + '.aln'

        reFilter = re.compile('Seq(\d*) *(.*)')
        #returns the 'n' of Seqn and the subsequent alignment string

        alnHandle = open(alnName, 'r+t')
        for line in alnHandle:
            reOutput = reFilter.findall(line)
            if len(reOutput) == 1:
                self.alignedseqs[int(reOutput[0][0])] += reOutput[0][1]


        alnHandle.close()

    def CleanUpFiles(self, FASTA = True, ALN = True, DND = True):
        """
        CleanUpFiles(fasta=True,aln=True,dnd=True)
            Cleans up the .fasta; .aln; and .dnd files that are left over
        after the alignment finishes.

        """
        if FASTA == True:
            try:
                os.remove(self.fastaname)
            except WindowsError:
                print 'Fasta Already Removed'
        if ALN:
            try:
                os.remove(self.fastaname[0:-6] + '.aln')
            except WindowsError:
                print 'ALN Already Removed'
        if DND:
            try:
                os.remove(self.fastaname[0:-6] + '.dnd')
            except WindowsError:
                print 'DND Already Removed'
        self.cleanedup = FASTA&ALN&DND
        
    def AlignSeqs(self, INPUT_SEQS = None, PARAM_PAIRS = None):
        """
        AlignSeqs(InputSeqs=None,paramPairs=None)
            A gateway function which takes a set of Input Sequences and
        Parameter Pairs and performs the Writing, Alignment, and Reading
        in one Call.

        """

        if INPUT_SEQS != None:
            self.CleanUpFiles(aln = False, dnd = False)
            self.fastawritten = False
            self.inputseqs = INPUT_SEQS

        if PARAM_PAIRS != None:
            self.parampairs = PARAM_PAIRS

        if not(self.fastawritten):
            self.WriteFASTA()
            
        self.DoAlignment()
        self.ReadALN()
        self.CleanUpFiles()

        return self.alignedseqs
        


