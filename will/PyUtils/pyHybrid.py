import subprocess
import os
import re
from Queue import *
import string
import threading
from HIVGenoTypeChecker import FastaDirIter
import time
import logging
import cPickle as pickle

DIREC = 'C:\\Documents and Settings\\Will\\My Documents\\HIVRNAi\\'
databaseName = 'human_predictions_4dwnload.txt'

logging.basicConfig(level=logging.DEBUG,format='[%(threadName)-10s] %(message)s %(asctime)s')


def ReadData(fileName):
    """
    ReadData
        Reads the miRNA database from the website: http://www.microrna.org/microrna/getDownloads.do

        miRNADict=ReadData(fileName)

            fileName        The full path to the database file.

            miRNADict       A dictionary where the KEY is a miRNA name and the value is the miRNA sequence
            

    """
    miRNAList = {}

    fileHandle = open(fileName,mode = 'r')

    #get rid of comment line
    thisLine = fileHandle.next()

    splitter = re.compile('\t')

    for thisLine in fileHandle:
        data = splitter.split(thisLine)
        miRNAList[data[3]] = string.upper(string.replace(data[5],'-',''))


    fileHandle.close()

    return miRNAList

def Calibrate(Seq,NUM_SAMPLES=5000, PRECOMPUTE_DICT = None,miRNAname = None,
              DUMP_QUEUE = None,CODING_FILE='C:\\RNAHybrid\\coding_HIV.txt',
              RNAHybridPath = 'C:\\RNAHybrid\\',PICKLE_DICT = None):
    """
    Calibrate
        Uses the RNAcalibrate function to determine the (delta,theta) of the extreme value
        distribution required for the RNAhybrid calculation.

    Calibrate(Seq, NUM_SAMPLES = 5000, PRECOMPUTE_DICT = None, miRNAname = None,
              DUMP_QUEUE = None,CODING_FILE='C:\\RNAHybrid\\coding_HIV.txt',
              RNAHybridPath = 'C:\\RNAHybrid\\',PICKLE_DICT = None)

              Seq           The miRNA sequence to calibrate

              NUM_SAMPLES   The number of random sequences to generate while determining
                            the extreme value distribution.

              PRECOMPUTE_DICT   A dictionary in which the KEY is a miRNA sequence and the
                                value is the (delta,theta) computed from a previous run.

              DUMP_QUEUE    A Queue object to dump the resulting output tuple

              miRNAName     The name of the miRNA being calibrated ... needed for DUMP_QUEUE

              CODING_FILE   The full path to a file which contains the dinucleotide
                            frequencies for the background sequence.

              RNAHybridPath The path to the directory where the RNAcalibrate.exe file is contained.

              PICKLE_DICT   The full path to a file in which to pickle a copy of the dict.

    RETURNS:

        tuple:
            (delta,theta) --- as STRINGS

        If DUMP_QUEUE != None
            The tuple is .put() onto the end of the queue
            (miRNAname, delta, theta) --- as STRINGS
    """

    #If a pre-computed dictionary is provided, check to see if it has the data we need
    try:
        outputList = ('junk','junk',PRECOMPUTE_DICT[miRNAname][0],PRECOMPUTE_DICT[miRNAname][1])
        logging.debug('Found Calib Data')
    except:
        logging.debug('Making Calib Data')
        command = RNAHybridPath + 'RNAcalibrate'

        command += ' -d ' + CODING_FILE
        command += ' -k ' + str(NUM_SAMPLES)

        command += ' ' + Seq

        SysCall = subprocess.Popen(command,shell=True,stdout = subprocess.PIPE)

        logging.debug('starting waiting')
        SysCall.wait()
        logging.debug('done waiting')

        output = SysCall.communicate()[0]
        outputList = re.split(' |\n',output)

    if DUMP_QUEUE == None:
        return (outputList[2],outputList[3])
    else:
        DUMP_QUEUE.put((miRNAname,Seq,outputList[2],outputList[3]))

    if PRECOMPUTE_DICT != None:
        PRECOMPUTE_DICT[miRNAname] = (outputList[2],outputList[3])

    if PICKLE_DICT != None:
        logging.debug('pickling Data')
        pickHandle = open(RNAHybridPath+'SavedCalib.pkl',mode='w')
        pickle.dump(PRECOMPUTE_DICT,pickHandle)
        pickHandle.close()
        
    
def HybridSeq(RNAiSeq,DeltaTheta,ChromSeq, miRNAname = None, ChomName = None, DUMP_QUEUE = None,RNAHybridPath = 'C:\\RNAHybrid\\'):
    """
    HybridSeq
        Makes a call to the RNAhybrid.exe to determine the likilyhood of the miRNA sequence binding to the Chomrosome
        Sequence.  It uses the data from Calibrate to determine the statistical significance.

        HybridSeq(RNAiSeq,DeltaTheta,ChromSeq, miRNAname = None, ChomName = None, DUMP_QUEUE = None,RNAHybridPath = 'C:\\RNAHybrid\\'):

        RNAiSeq         The sequence of the short miRNA.

        DeltaTheta      The output of Calibrate ... the shape of the extreme value distribution.

        ChromSeq        The long sequence of DNA or RNA to search for matches.

        miRNAname       The name of the miRNA being tested
        ChromName       The name of the chromosomal region being tested.
        
        DUMP_QUEUE      A queue to dump the resulting data into.

        RNAHybridPath   The path to the folder that contains RNAhybrid


    RETURNS

        the list of tuples:
            (position, energy, p-value)

        if DUMP_QUEUE != None
        The the following tuple is .put onto the Queue
            (ChomName, miRNAname, position, energy, p-value)


    """
    command = RNAHybridPath + 'RNAhybrid'

    command += ' -c'
    command += ' -d ' + DeltaTheta[0] + ',' + DeltaTheta[1]
    #command += ' -s ' + '3utr_human'
    command += ' -p 0.1' 
    command += ' ' + ChromSeq
    command += ' ' + RNAiSeq

    #logging.debug('Doing Hybrid')
    SysCall = subprocess.Popen(command,shell=True,stdout = subprocess.PIPE)
    SysCall.wait()
    

    output = SysCall.communicate()[0]
    
    finalOutput = []
    if len(output)>1:
        allLines = re.split('\n',output)
	#print len(allLines)
        for thisLine in allLines[0:-1]:
            #print thisLine
            finalOutput.append(re.match('(.*?):.*?command_line:[\-\.\d]*:([\.\-\d]*):([\.\-\d]*):([\.\-\d]*).*',thisLine).groups()[1:])

    if DUMP_QUEUE==None:
        return finalOutput
    else:
        if len(finalOutput)>0:
            for thisTup in finalOutput:
                DUMP_QUEUE.put((ChomName,miRNAname,thisTup[2],thisTup[0],thisTup[1]))
    

def LoadSeqs(LOAD_QUEUE,DUMP_QUEUE,AllCalibFlag,FinishedLoadFlag,fastaDir,oneDone):
    """
    LoadSeqs
        Uses the FastaDirIter to match each miRNA (and it DeltaTheta value) from the LOAD_QUEUE
        to the many sequence fragments returned by FastaDirIter.  It .put these onto the DUMP_QUEUE.

    LoadSeqs(LOAD_QUEUE,DUMP_QUEUE,AllCalibFlag,FinishedLoadFlag,fastaDir,oneDone)

        LOAD_QUEUE          A queue which is supplied by Calibrate where each element is a tuple:
                                    (miRNAname, delta, theta)

        DUMP_QUEUE          The place where each element is appended as a tuple:
                                    (ChromName, miRNAname, miRNAseq, ChromSeq, (delta, theta))

        AllCalibFlag        An event handle that is .set() when all calibrations have been completed
                            after which this function can exit.

        FinishedLoadFlag    An event handle which will be .set() when all calibrations and sequences
                            have been loaded into the DUMP_QUEUE

        fastaDir            The path to a directory with a collection of FASTA files.

        oneDone             An event handle which will be .set() when at least one element has been
                            loaded into the DUMP_QUEUE
                            

    """
    while True:
        try:
            thisCalib = LOAD_QUEUE.get_nowait()
        except:
            if AllCalibFlag.isSet():
                FinishedLoadFlag.set()
                break
            else:
                #print 'Didnt wait long enough'
                continue
        LOAD_QUEUE.task_done()
        logging.debug('Loading Seq')
        oneDone.set()
        for thisRec in FastaDirIter(fastaDir):
            
            DUMP_QUEUE.put((thisRec.id,thisCalib[0],thisCalib[1],thisRec.seq.tostring(),(thisCalib[2],thisCalib[3])))

        

        
        
def CalibWorker(LOAD_QUEUE,DumpingQueue,finishedEvent,oneDone,PreCompDict):
    """
    CalibWorker
        A worker function which facilitates the threading of Calibrate.

    CalibWorker(LOAD_QUEUE,DumpingQueue,finishedEvent,oneDone,PreCompDict)

        LOAD_QUEUE      A queue which contains tuples to be passed into Calibrate
                            (miRNAname, miRNAseq)

        DumpingQueue    A queue where Calibrate will dump its results.

        finishedEvent   An event handle which is .set() when the calibrating queue has been
                        completely emptied.

        oneDone         An event handle which is .set() when at least one element has been
                        loaded into the DumpingQueue

        PreCompDict     A dictionary which holds any pre-computed Calibrations


    """
    while True:
        try:
            thisMiRNA=LOAD_QUEUE.get_nowait()
        except:
           #print 'Finished Calibrating'
            finishedEvent.set()
            break
       #print 'Calibrating: ' + thisMiRNA[0]
        Calibrate(thisMiRNA[1],miRNAname=thisMiRNA[0],DUMP_QUEUE = DumpingQueue, PRECOMPUTE_DICT = PreCompDict)
        oneDone.set()
        LOAD_QUEUE.task_done()
        
def HybridWorker(LOAD_QUEUE,DumpingQueue,finishedLoad,finishedHybrid,oneDone):
    """
    HybridWorker
        A function which facilitates the threading of Hybrid.

    HybridWorker(LOAD_QUEUE,DumpingQueue,finishedLoad,finishedHybrid,oneDone)

        LOAD_QUEUE      A queue which has the tuples:
                            (ChromName,miRNAname,miRNAseq,ChromSeq,DeltaTheta)

        DumpingQueue    A queue to which Hybrid should dump is 'writable' data

        finishedLoad    An event handle which is .set() by LoadSeqs when all needed
                        sequences have been loaded into the queue.

        finishedHybrid  An event handle which is .set() when all of the hybridizations
                        have been completed.

        oneDone         An event handle which is .set() at least one Hyrbid has been
                        finished and loaded into the DumpingQueue

    """
    while True:
        try:
            thisSet = LOAD_QUEUE.get_nowait()
            #(ChromName,miRNAname,miRNAseq,ChromSeq,DeltaTheta)
        except:
            if finishedLoad.isSet():
               #print 'Finished Hybridizing'
                finishedHybrid.set()
                break
            else:
                #print 'Didnt wait long enough: HybridWorker'
                continue
       #print 'Hybing: ' + thisSet[0] + ' ' + thisSet[1]
        HybridSeq(thisSet[2],thisSet[4],thisSet[3],miRNAname=thisSet[1],ChomName=thisSet[0],DUMP_QUEUE=DumpingQueue)
        oneDone.set()
        LOAD_QUEUE.task_done()

def WritingWorker(LOAD_QUEUE,fHandle,finishedHybrid):
    """
    WritingWorker
        A function which writes the data from the LOAD_QUEUE into a file.

    WritingWorker(LOAD_QUEUE,fHandle,finishedHybrid)

        LOAD_QUEUE      A queue of tuples

        fHandle         The handle to an open file where the data is written too.
    
        finishedHybrid  An event handle which is .set() when all data has been processed
                        by the HybridWorker function.  This implies that when the LOAD_QUEUE
                        is empty then all data is complete.


    """
    while True:
        try:
            thisSet = LOAD_QUEUE.get_nowait()
        except:
            if finishedHybrid.isSet():
               #print 'Finished Writing'
                break
            else:
                #print 'Didt wait long enough: WritingWorker'
                continue
        
    #WriteData(fileName,thisSet)
        logging.debug('writting Data')
        fHandle.write(string.join(thisSet,sep='\t') + '\n')
        LOAD_QUEUE.task_done()


def DoAll(miRNADatabase,outputFile,NUM_CALIB = 1, NUM_HYBRID = 30):
    """
    DoAll
        A main function which initializes all of headwork and then begins the required threads.


    """
    miRNAdata = ReadData(miRNADatabase)

    try:
        pickHandle = open('C:\\RNAHybrid\\'+'SavedCalib.pkl',mode='r')
        preComputedCalibs = pickle.load(pickHandle)
        pickHandle.close()
        print 'Found SavedCalib File'
    except:
        print 'Could not Find SavedCalib File'
        preComputedCalibs={}
    
    NeedCalibQueue = Queue(-1)
    finishedCalibQueue = Queue(-1)
    NeedHybridQueue = Queue(500)
    NeedWriteQueue = Queue(500)

    counter=0
    for thisKey in miRNAdata:
        NeedCalibQueue.put_nowait((thisKey,miRNAdata[thisKey]))


    finishedCalib = threading.Event()
    oneCalib = threading.Event()
    finishedLoading = threading.Event()
    oneLoad = threading.Event()
    finishedHybriding = threading.Event()
    oneHybrid = threading.Event()
    finishedWriting = threading.Event()

    #calibThreadVec = []
    for i in range(NUM_CALIB):
        calibThread = threading.Thread(name='calibThread',target=CalibWorker,args=(NeedCalibQueue,finishedCalibQueue,finishedCalib,oneCalib,preComputedCalibs))
    #calibThreadVec.append(calibThread)
        calibThread.start()
    #CalibWorker(NeedCalibQueue,finishedCalibQueue,finishedCalib)    

    oneCalib.wait()
    LoadSeqThread = threading.Thread(name='LoadSeqThread',target=LoadSeqs,args=(finishedCalibQueue,NeedHybridQueue,finishedCalib,finishedLoading,'C:\\RNAHybrid\\seqs\\',oneLoad))
    LoadSeqThread.start()
    #LoadSeqs(finishedCalibQueue,NeedHybridQueue,finishedCalib,finishedLoading,'C:\\RNAHybrid\\seqs\\',oneLoad)

    oneLoad.wait()
    HybridThreadVec = []
    for i in range(NUM_HYBRID):
        HybridThread = threading.Thread(name='HybridThread' + str(i),target=HybridWorker,args=(NeedHybridQueue,NeedWriteQueue,finishedLoading,finishedHybriding,oneHybrid))
        HybridThreadVec.append(HybridThread)
        HybridThread.start()
        
    #HybridWorker(NeedHybridQueue,NeedWriteQueue,finishedLoading,finishedHybriding)

    oneHybrid.wait()
    outputHandle = open(outputFile,'w+')
    WrittingThread = threading.Thread(name='WrittingThread',target=WritingWorker,args=(NeedWriteQueue,outputHandle,finishedHybriding))
    WrittingThread.start()
    #WritingWorker(NeedWriteQueue,outputHandle,finishedHybriding)

    WrittingThread.join()
    outputHandle.close()
    




if __name__=='__main__':
    DoAll(DIREC+databaseName,'C:\\RNAHybrid\\newoutput.txt')


