"""
pyHybrid
    A set of functions for interfacing with the RNAhybrid program.
"""

import subprocess
import re
import Queue
import string
import threading
from HIVGenoTypeChecker import FastaDirIter
import logging
import cPickle as pickle

DIREC = 'C:\\Documents and Settings\\Will\\My Documents\\HIVRNAi\\'
DATABASENAME = 'human_predictions_4dwnload.txt'

logging.basicConfig(level = logging.DEBUG,
                    format = '[%(threadName)-10s] %(message)s %(asctime)s')


def ReadData(FILENAME):
    """
    ReadData
        Reads the miRNA database from the website:
            http://www.microrna.org/microrna/getDownloads.do

        miRNADict=ReadData(fileName)

            fileName        The full path to the database file.

            miRNADict       A dictionary where the KEY is a miRNA name and the
                            value is the miRNA sequence
            

    """
    mi_rna_list = {}

    file_handle = open(FILENAME, mode = 'r')

    #get rid of comment line
    this_line = file_handle.next()

    splitter = re.compile('\t')

    for this_line in file_handle:
        data = splitter.split(this_line)
        mi_rna_list[data[3]] = string.upper(string.replace(data[5], '-', ''))


    file_handle.close()

    return mi_rna_list

def Calibrate(SEQ, NUM_SAMPLES = 5000, PRECOMPUTE_DICT = None,
              MI_RNA_NAME = None, DUMP_QUEUE = None,
              CODING_FILE = 'C:\\RNAHybrid\\coding_HIV.txt',
              RNAHybridPath = 'C:\\RNAHybrid\\', PICKLE_DICT = None):
    """
    Calibrate
        Uses the RNAcalibrate function to determine the (delta,theta) of the
        extreme value distribution required for the RNAhybrid calculation.

    Calibrate(SEQ, NUM_SAMPLES = 5000, PRECOMPUTE_DICT = None,
              MI_RNA_NAME = None,
              DUMP_QUEUE = None,CODING_FILE='C:\\RNAHybrid\\coding_HIV.txt',
              RNAHybridPath = 'C:\\RNAHybrid\\',PICKLE_DICT = None)

              SEQ           The miRNA sequence to calibrate

              NUM_SAMPLES   The number of random sequences to generate
                            while determining the extreme value distribution.

              PRECOMPUTE_DICT   A dictionary in which the KEY is a miRNA
                                sequence and the value is the (delta,theta)
                                computed from a previous run.

              DUMP_QUEUE    A Queue object to dump the resulting output tuple

              miRNAName     The name of the miRNA being
                            calibrated ... needed for DUMP_QUEUE

              CODING_FILE   The full path to a file which contains the
                            dinucleotide frequencies for the background
                            sequence.

              RNAHybridPath The path to the directory where the
                            RNAcalibrate.exe file is contained.

              PICKLE_DICT   The full path to a file in which to pickle a copy
                            of the dict.

    RETURNS:

        tuple:
            (delta,theta) --- as STRINGS

        If DUMP_QUEUE != None
            The tuple is .put() onto the end of the queue
            (MI_RNA_NAME, delta, theta) --- as STRINGS
    """

    #If a pre-computed dictionary is provided,
    #check to see if it has the data we need
    try:
        outputList = ('junk', 'junk', PRECOMPUTE_DICT[MI_RNA_NAME][0],
                      PRECOMPUTE_DICT[MI_RNA_NAME][1])
        logging.debug('Found Calib Data')
    except KeyError:
        logging.debug('Making Calib Data')
        command = RNAHybridPath + 'RNAcalibrate'

        command += ' -d ' + CODING_FILE
        command += ' -k ' + str(NUM_SAMPLES)

        command += ' ' + SEQ

        sys_call = subprocess.Popen(command, shell = True,
                                   stdout = subprocess.PIPE)

        logging.debug('starting waiting')
        sys_call.wait()
        logging.debug('done waiting')

        output = sys_call.communicate()[0]
        outputList = re.split(' |\n', output)

    if DUMP_QUEUE == None:
        return (outputList[2], outputList[3])
    else:
        DUMP_QUEUE.put((MI_RNA_NAME, SEQ, outputList[2], outputList[3]))

    if PRECOMPUTE_DICT != None:
        PRECOMPUTE_DICT[MI_RNA_NAME] = (outputList[2], outputList[3])

    if PICKLE_DICT != None:
        logging.debug('pickling Data')
        pickHandle = open(RNAHybridPath + 'SavedCalib.pkl', mode = 'w')
        pickle.dump(PRECOMPUTE_DICT, pickHandle)
        pickHandle.close()
        
    
def HybridSeq(RNA_SEQ, DELTA_THETA, CHROM_SEQ, MI_RNA_NAME = None,
              CHROM_NAME = None, DUMP_QUEUE = None,
              RNA_HYBRID_PATH = 'C:\\RNAHybrid\\'):
    """
    HybridSeq
        Makes a call to the RNAhybrid.exe to determine the likilyhood of the
        miRNA sequence binding to the Chomrosome Sequence.  It uses the data
        from Calibrate to determine the statistical significance.

        HybridSeq(RNA_SEQ,DELTA_THETA,CHROM_SEQ, MI_RNA_NAME = None,
              CHROM_NAME = None, DUMP_QUEUE = None,
              RNA_HYBRID_PATH = 'C:\\RNAHybrid\\'):

        RNA_SEQ         The sequence of the short miRNA.

        DELTA_THETA      The output of Calibrate ... the shape of the extreme
                         value distribution.

        CHROM_SEQ        The long sequence of DNA or RNA to search for matches.

        MI_RNA_NAME       The name of the miRNA being tested
        ChromName       The name of the chromosomal region being tested.
        
        DUMP_QUEUE      A queue to dump the resulting data into.

        RNA_HYBRID_PATH   The path to the folder that contains RNAhybrid


    RETURNS

        the list of tuples:
            (position, energy, p-value)

        if DUMP_QUEUE != None
        The the following tuple is .put onto the Queue
            (ChomName, MI_RNA_NAME, position, energy, p-value)


    """
    command = RNA_HYBRID_PATH + 'RNAhybrid'

    command += ' -c'
    command += ' -d ' + DELTA_THETA[0] + ',' + DELTA_THETA[1]
    #command += ' -s ' + '3utr_human'
    command += ' -p 0.1' 
    command += ' ' + CHROM_SEQ
    command += ' ' + RNA_SEQ

    #logging.debug('Doing Hybrid')
    sys_call = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE)
    sys_call.wait()
    

    output = sys_call.communicate()[0]

    output_reg_exp = '(.*?):.*?command_line:[\-\.\d]*:([\.\-\d]*):'
    output_reg_exp += '([\.\-\d]*):([\.\-\d]*).*'
    
    final_output = []
    if len(output) > 1:
        all_lines = re.split('\n', output)
	#print len(allLines)
        for this_line in all_lines[0:-1]:
            #print this_line
            final_output.append(re.match(output_reg_exp,
                                         this_line).groups()[1:])

    if DUMP_QUEUE == None:
        return final_output
    else:
        if len(final_output) > 0:
            for thistup in final_output:
                DUMP_QUEUE.put((CHROM_NAME, MI_RNA_NAME, thistup[2],
                                thistup[0], thistup[1]))
    

def LoadSeqs(LOAD_QUEUE, DUMP_QUEUE, ALL_CALIB_FLAG, FINISHED_LOAD_FLAG,
             FASTA_DIR, ONE_DONE):
    """
    LoadSeqs
        Uses the FastaDirIter to match each miRNA (and it DELTA_THETA value)
        from the LOAD_QUEUE to the many sequence fragments returned by
        FastaDirIter.  It .put these onto the DUMP_QUEUE.

    LoadSeqs(LOAD_QUEUE,DUMP_QUEUE,ALL_CALIB_FLAG,FINISHED_LOAD_FLAG,
             FASTA_DIR,ONE_DONE)

        LOAD_QUEUE          A queue which is supplied by Calibrate where each
                            element is a tuple:
                                    (MI_RNA_NAME, delta, theta)

        DUMP_QUEUE          The place where each element is appended as a tuple:
                                    (ChromName, MI_RNA_NAME, miRNAseq,
                                    CHROM_SEQ, (delta, theta))

        AllCalibFlag        An event handle that is .set() when all calibrations
                            have been completed after which this function
                            can exit.

        FinishedLoadFlag    An event handle which will be .set() when all
                            calibrations and sequences have been loaded into
                            the DUMP_QUEUE

        fastaDir            The path to a directory with a collection of
                            FASTA files.

        oneDone             An event handle which will be .set() when at least
                            one element has been loaded into the DUMP_QUEUE
                            

    """
    while True:
        try:
            this_calib = LOAD_QUEUE.get_nowait()
        except  :
            if ALL_CALIB_FLAG.isSet():
                FINISHED_LOAD_FLAG.set()
                break
            else:
                #print 'Didnt wait long enough'
                continue
        LOAD_QUEUE.task_done()
        logging.debug('Loading SEQ')
        ONE_DONE.set()
        for this_rec in FastaDirIter(FASTA_DIR):
            
            DUMP_QUEUE.put((this_rec.id, this_calib[0], this_calib[1],
                            this_rec.seq.tostring(),
                            (this_calib[2], this_calib[3])))

        

        
        
def CalibWorker(LOAD_QUEUE, DUMPING_QUEUE, FINISHED_EVENT,
                ONE_DONE, PRE_COMP_DICT):
    """
    CalibWorker
        A worker function which facilitates the threading of Calibrate.

    CalibWorker(LOAD_QUEUE,DUMPING_QUEUE,FINISHED_EVENT,oneDone,PRE_COMP_DICT)

        LOAD_QUEUE      A queue which contains tuples to be passed into
                        Calibrate:
                            (MI_RNA_NAME, miRNAseq)

        DUMPING_QUEUE    A queue where Calibrate will dump its results.

        FINISHED_EVENT   An event handle which is .set() when the calibrating
                         queue has been completely emptied.

        oneDone         An event handle which is .set() when at least one
                        element has been loaded into the DUMPING_QUEUE

        PRE_COMP_DICT     A dictionary which holds any pre-computed Calibrations


    """
    while True:
        try:
            this_mirna = LOAD_QUEUE.get_nowait()
        except  :
           #print 'Finished Calibrating'
            FINISHED_EVENT.set()
            break
       #print 'Calibrating: ' + this_mirna[0]
        Calibrate(this_mirna[1], MI_RNA_NAME = this_mirna[0],
                  DUMP_QUEUE = DUMPING_QUEUE, PRECOMPUTE_DICT = PRE_COMP_DICT)
        ONE_DONE.set()
        LOAD_QUEUE.task_done()
        
def HybridWorker(LOAD_QUEUE, DUMPING_QUEUE, FINISHED_LOAD,
                 FINISHED_HYBRID, ONE_DONE):
    """
    HybridWorker
        A function which facilitates the threading of Hybrid.

    HybridWorker(LOAD_QUEUE,DUMPING_QUEUE,FINISHED_LOAD,
                 FINISHED_HYBRID,ONE_DONE):

        LOAD_QUEUE      A queue which has the tuples:
                        (ChromName,MI_RNA_NAME,miRNAseq,CHROM_SEQ,DELTA_THETA)

        DUMPING_QUEUE    A queue to which Hybrid should dump is 'writable' data

        FINISHED_LOAD    An event handle which is .set() by LoadSeqs when all
                        needed sequences have been loaded into the queue.

        FINISHED_HYBRID  An event handle which is .set() when all of the
                         hybridizations have been completed.

        oneDone         An event handle which is .set() at least one Hyrbid has
                        been finished and loaded into the DUMPING_QUEUE

    """
    while True:
        try:
            this_set = LOAD_QUEUE.get_nowait()
            #(ChromName,MI_RNA_NAME,miRNAseq,CHROM_SEQ,DELTA_THETA)
        except  :
            if FINISHED_LOAD.isSet():
               #print 'Finished Hybridizing'
                FINISHED_HYBRID.set()
                break
            else:
                #print 'Didnt wait long enough: HybridWorker'
                continue
       #print 'Hybing: ' + this_set[0] + ' ' + this_set[1]
        HybridSeq(this_set[2], this_set[4], this_set[3],
                  MI_RNA_NAME = this_set[1], CHROM_NAME = this_set[0],
                  DUMP_QUEUE = DUMPING_QUEUE)
        
        ONE_DONE.set()
        LOAD_QUEUE.task_done()

def WritingWorker(LOAD_QUEUE, FILE_HANDLE, FINISHED_HYBRID):
    """
    WritingWorker
        A function which writes the data from the LOAD_QUEUE into a file.

    WritingWorker(LOAD_QUEUE,FILE_HANDLE,FINISHED_HYBRID)

        LOAD_QUEUE      A queue of tuples

        FILE_HANDLE         The handle to an open file where the data is
                            written too.
    
        FINISHED_HYBRID  An event handle which is .set() when all data has
                         been processed by the HybridWorker function.  This
                         implies that when the LOAD_QUEUE is empty then all
                         data is complete.


    """
    while True:
        try:
            this_set = LOAD_QUEUE.get_nowait()
        except  :
            if FINISHED_HYBRID.isSet():
               #print 'Finished Writing'
                break
            else:
                #print 'Didt wait long enough: WritingWorker'
                continue
        
    #WriteData(fileName,this_set)
        logging.debug('writting Data')
        FILE_HANDLE.write(string.join(this_set, sep = '\t') + '\n')
        LOAD_QUEUE.task_done()


def DoAll(MI_RNA_DATABASE, OUTPUT_FILE, NUM_CALIB = 1, NUM_HYBRID = 30):
    """
    DoAll
        A main function which initializes all of headwork and then begins the
        required threads.


    """
    mi_rna_data = ReadData(MI_RNA_DATABASE)

    try:
        pickHandle = open('C:\\RNAHybrid\\' + 'SavedCalib.pkl', mode = 'r')
        preComputedCalibs = pickle.load(pickHandle)
        pickHandle.close()
        print 'Found SavedCalib File'
    except WindowsError:
        print 'Could not Find SavedCalib File'
        preComputedCalibs={}
    
    need_calib_queue = Queue.Queue(-1)
    finished_calib_queue = Queue.Queue(-1)
    need_hybrid_queue = Queue.Queue(500)
    need_write_queue = Queue.Queue(500)

    for this_hey in mi_rna_data:
        need_calib_queue.put_nowait((this_hey, mi_rna_data[this_hey]))


    finished_calib = threading.Event()
    one_calib = threading.Event()
    finished_loading = threading.Event()
    one_load = threading.Event()
    finished_hybriding = threading.Event()
    one_hybrid = threading.Event()


    for i in range(NUM_CALIB):
        cailib_thread = threading.Thread(name = 'cailib_thread',
                                         target = CalibWorker,
                                         args = (need_calib_queue,
                                                 finished_calib_queue,
                                                 finished_calib, one_calib,
                                                 preComputedCalibs))

        cailib_thread.start()


    one_calib.wait()
    load_seq_thread = threading.Thread(name = 'load_seq_thread',
                                       target = LoadSeqs,
                                       args = (finished_calib_queue,
                                               need_hybrid_queue,
                                               finished_calib, finished_loading,
                                               'C:\\RNAHybrid\\seqs\\',
                                               one_load))
    load_seq_thread.start()


    one_load.wait()
    hybrid_thread_vec = []
    for i in range(NUM_HYBRID):
        hybrid_thread = threading.Thread(name = 'hybrid_thread' + str(i),
                                         target = HybridWorker,
                                         args = (need_hybrid_queue,
                                                 need_write_queue,
                                                 finished_loading,
                                                 finished_hybriding,
                                                 one_hybrid))
        
        hybrid_thread_vec.append(hybrid_thread)
        hybrid_thread.start()
        
    one_hybrid.wait()
    output_handle = open(OUTPUT_FILE,'w+')
    writting_thread = threading.Thread(name = 'writting_thread',
                                      target = WritingWorker,
                                      args = (need_write_queue, output_handle,
                                              finished_hybriding))
    
    writting_thread.start()

    writting_thread.join()
    output_handle.close()
    




if __name__ == '__main__':
    DoAll(DIREC + DATABASENAME, 'C:\\RNAHybrid\\newoutput.txt')


