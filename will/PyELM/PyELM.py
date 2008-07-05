import os
import numpy
import re
from pyQueueUtils import *
from copy import *
import threading
import time
import bisect
import Queue
import types
import itertools


def ParseELMs(SEQ, ELM_DICT):
    """
    ParseELMs
        Searches for the regular expressions of ELM_DICT within the
        SEQUENCES provided.

    RESULT_DICT = ELMMatchSeqs(SEQ,ELM_DICT)

    SEQ             The sequence to be parsed.
    ELM_DICT        A dict keyed by ELM names and has a value of re.compile


    RESULT_DICT     A map keyed by ELM names with a numpy array of matched
                    positions

    """
    
    seq_elm_dict = {}

    for this_elm in ELM_DICT:
        this_index = []
        #each search only finds one match
        #so repeat the search until no more are found
        spot = ELM_DICT[this_elm].search(SEQ)
        while spot != None:
            this_index.append(spot.start())
            spot = ELM_DICT[this_elm].search(SEQ, spot.start() + 1)
        seq_elm_dict[this_elm] = numpy.array(this_index)
    
    return seq_elm_dict

class SeqObject:
    def __init__(self, SEQ, ELM_DICT = None):
        self.sequence = SEQ

        if ELM_DICT != None:
            self.elm_positions = ParseELMs(self.sequence, ELM_DICT)


    def GetResult(self, WANTED_ELM, RANGE, OUTPUT_TYPE = 'count'):
        """
        GetResult
            Returns the results of of the ELM parsing in various formats based
            on the input sets.

        OUTPUT = GetResult(WANTED_ELM, RANGE, OUTPUT_TYPE)

        RANGE       OUTPUT_TYPE         OUTPUT
        'None'      'count'             The total number of WANTED_ELM matches
                                        along the sequence

        [B E]       'count'             The number of WANTED_ELM matches between
                                        B and E

        [B E]       'array'             A boolean array where each spot
                                        indicates a presence/absence call.
        """

        if RANGE == None:
            RANGE = (0, len(self.sequence))

        this_array = self.elm_positions[WANTED_ELM]

        if OUTPUT_TYPE == 'count':
            return numpy.sum(numpy.logical_and(this_array >= RANGE[0],
                                               this_array <= RANGE[1]))
        
        if OUTPUT_TYPE == 'array':
            valid_spots = numpy.logical_and(this_array >= RANGE[0],
                                            this_array <= RANGE[1])
            modified_spots = this_array[valid_spots] - RANGE[0]
            output_array = numpy.zeros((RANGE[1] - RANGE[0], 1))
            output_array[modified_spots] = 1
            return output_array

    def AddELM(self, ELM_NAME, REG_EXPR):
        """
        AddELM
            Adds new ELMs to the parsed database.

            AddELM(ELM_NAME, REG_EXPR)

        """

        temp_dict = {}
        temp_dict[ELM_NAME] = REG_EXPR

        self.elm_positions.update(ParseELMs(self.sequence,temp_dict))
        


class PyELM:

    def __init__(self):

        self.elm_dict = {}
        self.elm_order = []
        self.sequences = []
        self.elm_re_dict = {}
        self.max_size = 0

        self.__found_correct = False

        self.precalc_elm_bin_dict = {}

        self.__singlebin_dict = {}
        self.__all_bin_dict = {}
        self.__depth_checked = []
        self.__working_queue = None

        self.__global_counter = 0
        self.max_breadth = 4
        self.max_depth = 100
        self.max_width = 5
        self.max_iter = 1000
        self.calc_time_out = 300
        self.force_bins = [858, 1359, 1483, 1676, 2570, 3434, 3513]

    def ELMParser(self, DIRECTORY = None):
        """
        Parser
            Parses a directory of HTML files to extract the ELM names and
            Regular Expressions and instanciates it into SELF.

            Parser(DIRECTORY)

            DIRECTORY       The directory containing the HTML files downloaded
                            from the ELM database.  If left empty then
            C:\Documents and Settings\Will\My Documents\PyELM\ELM_RAW_DOWNLOAD\

        """
        if DIRECTORY == None:
            DIRECTORY = "C:\\Documents and Settings\\Will\\My Documents\\ELM_Motif_finder\\ELM_RAW_DOWNLOAD\\"

        file_list = os.listdir(DIRECTORY)
        all_reg_exps = {}
        
        for i in xrange(1, len(file_list) - 1):
            
            this_file = open(DIRECTORY + file_list[i], 'r')
            elm_name = file_list[i].rpartition('.')
            self.elm_order.append(elm_name[0])
            current_line = ''
            
            while current_line.find('Pattern:') == -1:
                current_line = this_file.readline()
            
            current_line = this_file.readline()
            
            all_reg_exps[elm_name[0]] = current_line[current_line.find('>') +
                                                     1:current_line.find('/') - 2]
            this_file.close()
        
        self.elm_dict = all_reg_exps

    def AddELM(self, NAME, REG_EXP):
        """ 
        AddELM
            Loads a single ELM RegExp into the database.

        AddELM(Name,RegExp)

        """

        self.elm_order.append(NAME)
        self.elm_dict[NAME] = REG_EXP
        self.elm_re_dict[NAME] = re.compile(REG_EXP)

        for this_seq in self.sequences:
            this_seq.AddELM(NAME, self.elm_re_dict[NAME])
        
    def LoadSeqs(self, SEQUENCES):
        """
        LoadSeqs
            Loads a set of sequences into this instance for further analysis.

            LoadSeqs(SEQUENCES)

            SEQUENCES       A LIST of characters representing Amino Acids


        """

        
        for this_seq in SEQUENCES:
            self.sequences.append(SeqObject(this_seq, ELM_DICT = self.elm_re_dict))

    def GetSequence(self, SEQ_IND, RANGE = None):
        """
        GetSequence
            Returns the sequence held by the object.

        GetSequence(RANGE = None)
        """
        if RANGE == None:
            return self.sequences[SEQ_IND].sequence
        else:
            return self.sequences[SEQ_IND].sequence[RANGE[0]:RANGE[1]]

    def SimpleELMCall(self, SEQ_IND, WANTED_ELM,
                      WANTED_RANGE = None, OUTPUT_TYPE = 'bool'):
        """
        SimpleELMCall
            Returns a simple Presence/Absence Call for the ELM.

        PACall = SimpleELMCall(SEQ_IND,WANTED_ELM)

        SEQ_IND     The index of the desired sequence
        WANTED_ELM  The desired ELM

        PACall      Returns 1 if the ELM is present in the sequence and
                    0 otherwise.

        """
        if OUTPUT_TYPE == 'bool':
            return self.sequences[SEQ_IND].GetResult(WANTED_ELM, WANTED_RANGE) > 0
        elif OUTPUT_TYPE == 'count':
            return self.sequences[SEQ_IND].GetResult(WANTED_ELM, WANTED_RANGE)

    def MultiELMCall(self, SEQ_IND, WANTED_ELM):
        """
        MultiELMCall
            Returns the binned Presence/Absence Call for the ELM.

        PACall = MultiELMCall(SEQ_IND,WANTED_ELM)

        SEQ_IND     The index of the desired sequence
        WANTED_ELM  The desired ELM

        PACall      Returns 1 if the ELM is present in the sequence and
                    0 otherwise.

        """
        if WANTED_ELM not in self.precalc_elm_bin_dict:
            raise KeyError
            
        thisCall = []
        allBins = self.precalc_elm_bin_dict[WANTED_ELM]
        for i in xrange(len(allBins) - 1):
            if self.SimpleELMCall(SEQ_IND, WANTED_ELM,
                                  WANTED_RANGE = (allBins[i], allBins[i+1])):
                thisCall.append(1)
            else:
                thisCall.append(0)
        return thisCall

    def GetIter(self, ITER_TYPE, ITER_RETURN, ITER_RANGE):
        """
        SetIterMethod
            Sets the iterator method to the desired type and output.

        SetIterMethod(ITER_TYPE,ITER_RETURN)

        Potential Combinations:
            ITER_TYPE       ITER_RETURN     ITER_RANGE
        Method 1
            'ELM'           'Name'          None
            Iterates over the names in ELMDict
        Method 2
            'ELM'           [SeqNum]        None
            Iterates over the Presence/Absence calls of each ELM over the
            SeqNum'th Sequence
        Method 3    
            'SEQ'           'Sequence'      [RANGE]
            Iterates over each Sequence and returns the sequence from RANGE[0]
            to RANGE[1].  If RANGE == None then the entire sequence is returned
        Method 4    
            'SEQ'           'ELMname'       None
            Iterates over each sequence and returns P/A call of the ELMname
        Method 5    
            'SEQ'           'ELMname'       [RANGE]
            Iterates over each sequence and returns P/A call of ELMs found
            over desired range.
        Method 6
            'COUNT'         'ELMname'       None
            Iterates over each sequence and returns count of the ELMs in the
            sequence.
        Method 7
            'COUNT'         'ELMname'       [RANGE]
            Iterates over each sequence and returns count of the ELMs over the
            RANGE.

            

        """
        #Method 1
        if (ITER_TYPE == 'ELM') & (ITER_RETURN == 'Name'):
            return self.elm_order.__iter__()

        #Method 2
        if (ITER_TYPE == 'ELM') & (types.IntType == type(ITER_RETURN)):
            return itertools.imap(self.SimpleELMCall,
                                  itertools.repeat(ITER_RETURN),
                                  self.elm_order.__iter__())

        #Method 3
        if (ITER_TYPE == 'SEQ'):
            if (ITER_RETURN == 'Sequence'):
                return itertools.imap(self.GetSequence, itertools.count(),
                                      itertools.repeat(ITER_RANGE,
                                                       len(self.sequences)))
        #Method 4 & 5
            if (self.elm_dict.has_key(ITER_RETURN)):
                return itertools.imap(self.SimpleELMCall,
                                      itertools.count(),
                                      itertools.repeat(ITER_RETURN,
                                                       len(self.sequences)),
                                      itertools.repeat(ITER_RANGE))
        #Method 6 & 7
        if (ITER_TYPE == 'COUNT'):
            if (self.elm_dict.has_key(ITER_RETURN)):
                return itertools.imap(self.SimpleELMCall,
                                      itertools.count(),
                                      itertools.repeat(ITER_RETURN,
                                                       len(self.sequences)),
                                      itertools.repeat(ITER_RANGE),
                                      itertools.repeat('count'))
            

        raise KeyError

    def __QueueWorker(self, WANTED_ELM):
        """
        A helper function to manage the Priority Queue and Threads
        """
        while True:
            self.__global_counter += 1 
            try:
                #print 'Worker qsize: ', self.__working_queue.qsize()
                thisItem = self.__working_queue.get(True, 30)
                
            except Queue.Empty:
                #print 'queue empty'
                break
            
            if self.__found_correct.isSet():
                #print 'Detected Bailing, dumping stuff out of queue'
                self.__working_queue.task_done()
                continue

            self.__MakeBreadth(WANTED_ELM, thisItem[0],
                               thisItem[1], thisItem[2])
            self.__working_queue.task_done()
                

    def __MakeBreadth(self, WANTED_ELM, ALL_BINS, B_COUNTER, D_COUNTER):
        """
        __MakeBreadth

        Performs a breadth-first search.  All nodes are added to the Queue and
        processes in order.


        """
        
        if D_COUNTER > self.max_depth:
            return

        if (B_COUNTER > self.max_breadth) & (ALL_BINS in
                                                  self.__depth_checked):
            return
            
        local_q = []
        local_vals = []
            
        for k in xrange(len(ALL_BINS) - 1):
            bin_set = deepcopy(ALL_BINS)
            new_edge = bin_set[k] + int((bin_set[k + 1] - bin_set[k]) / 2)
            bin_set.insert(k + 1, new_edge)
            local_q.append(bin_set)
            local_vals.append(self.__EvaluateAllBins(WANTED_ELM, bin_set))

        for k in xrange(1, len(ALL_BINS) - 1):
            bin_set = deepcopy(ALL_BINS)
            bin_set.pop(k)
            local_q.append(bin_set)
            local_vals.append(self.__EvaluateAllBins(WANTED_ELM, bin_set))

        for k in xrange(1, len(ALL_BINS) - 1):
            bin_set = deepcopy(ALL_BINS)
            bin_move = int((bin_set[k+1]-bin_set[k])/2)
            local_q.append(bin_set)
            local_vals.append(self.__EvaluateAllBins(WANTED_ELM, bin_set))

        for k in xrange(1, len(ALL_BINS) - 1):
            bin_set = deepcopy(ALL_BINS)
            bin_move = int((bin_set[k] - bin_set[k - 1]) / 2)
            bin_set[k] = bin_set[k] - bin_move
            local_q.append(bin_set)
            local_vals.append(self.__EvaluateAllBins(WANTED_ELM, bin_set))

        best_inds = numpy.argsort(local_vals)


        if B_COUNTER < self.max_breadth:
            for k in xrange(min(self.max_width, len(best_inds)) - 1, 0, -1):
                self.__working_queue.put(((local_q[best_inds[k]],
                                           B_COUNTER + 1, D_COUNTER),
                                          local_vals[best_inds[k]]))
        else:
            if (local_vals[best_inds[0]] < self.__EvaluateAllBins(WANTED_ELM,
                                                                 ALL_BINS)):

                self.__depth_checked.append(ALL_BINS)
                self.__working_queue.put(((local_q.pop(best_inds[0]),
                                           B_COUNTER + 1, D_COUNTER+ 1),
                                          local_vals[best_inds[0]]))


    

    def __EvaluateBin(self, WANTED_ELM, BIN_EDGES):
        """
        Evaluates the values of the single bin provided.
        """
        bin_val = [0, 0, 0]
        for current_bin in self.GetIter('COUNT', WANTED_ELM,
                                        (BIN_EDGES[0], BIN_EDGES[1])):
            if current_bin == 0:
                bin_val[0] += 1
            elif current_bin == 1:
                bin_val[1] += 1
            else:
                bin_val[2] += 1
        return (bin_val[0], bin_val[1], bin_val[2] )


    def __EvaluateAllBins(self, WANTED_ELM, THESE_BINS):
        """
        Evaluates the values of the bins provided.  Uses a dictionary to avoid
        re-calculating values wherever possible
        """
        
        if len(self.force_bins) > 0:
            all_bins = deepcopy(THESE_BINS)
            for forced in self.force_bins:
                bisect.insort(all_bins,forced)
        else:
            all_bins = THESE_BINS
        
        if not(tuple(all_bins) in self.__all_bin_dict ):
            total_empty = 0
            total_correct = 0
            total_wrong = 0
            for i in xrange(len(all_bins) - 1):
                if (all_bins[i], all_bins[i + 1]) not in self.__singlebin_dict:
                    self.__singlebin_dict[(all_bins[i], all_bins[i + 1])] = self.__EvaluateBin(WANTED_ELM, (all_bins[i], all_bins[i + 1]))
                
                total_empty += self.__singlebin_dict[(all_bins[i],
                                                     all_bins[i + 1])][0]
                total_correct += self.__singlebin_dict[(all_bins[i],
                                                       all_bins[i + 1])][1]
                total_wrong += self.__singlebin_dict[(all_bins[i],
                                                     all_bins[i + 1])][2]
            funVal = (.25*total_empty + total_wrong) / (total_correct + 1)
            self.__all_bin_dict[tuple(all_bins)] = (total_empty, total_correct,
                                                    total_wrong, funVal)
            if (self.__all_bin_dict[tuple(all_bins)][2] == 0) & (self.__all_bin_dict[tuple(all_bins)][0] == 0):

                self.__found_correct.set()
        return self.__all_bin_dict[tuple(all_bins)][3]

    def GetAllBinError(self, WANTED_ELM, ALL_BINS):
        """
        Evaluates the values of the bins provided.
        """

        this_list = [];
                
        for i in xrange(len(ALL_BINS) - 1):
            this_bin=self.__EvaluateBin(WANTED_ELM,
                                       (ALL_BINS[i], ALL_BINS[i + 1]))
                        
            empty = this_bin[0]
            correct = this_bin[1]
            wrong = this_bin[2]

            this_list.append((empty, correct, wrong))    

        return this_list

    def __BailingThread(self):
        """
        __BailingThread
            Bails out of the current calculation.
        """
        print 'Exceded timeLimit, bailing on further exploration'
        self.__found_correct.set()
        return

    def CalculateBins(self, WANTED_ELM):
        """
        CalculateBins
            Uses a dynamic programing algorithm to find the optimal binning solution
        """
        num_thread = 4

        self.__singlebin_dict = {}
        self.__all_bin_dict = {}
        self.__depth_checked = []
        self.__found_correct = threading.Event()
        self.__global_counter = 0
        
        self.__working_queue = PriorityQueue(-1)
        these_bins = [0, int(self.max_size / 2), self.max_size]

        self.__working_queue.put(((these_bins, 0, 0), 0))

        worker_threads = []
        for i in range(num_thread-1):
            print 'Starting Thread: ', i
            t = threading.Thread(target = self.__QueueWorker,
                                 args = (WANTED_ELM,))
            worker_threads.append(t)
            time.sleep(15)
            t.start()

            
        print 'Starting Bailing Thread'
        stopping_thread = threading.Timer(self.calc_time_out,
                                         self.__BailingThread)
        stopping_thread.start()

        self.__working_queue.join()
        stopping_thread.cancel()

        
        for thisThread in worker_threads:
            print 'Joining Thread'
            thisThread.join()
        
        #print 'Calculating Values'
        
            
            

        #self.__QueueWorker(WANTED_ELM)

        currentMin = 50000
        currentMinInd = 0
        for i in self.__all_bin_dict:
            if self.__all_bin_dict[i][3] < currentMin:
                currentMinInd = i
                currentMin = self.__all_bin_dict[currentMinInd][3]

        self.precalc_elm_bin_dict[WANTED_ELM] = currentMinInd
        self.__working_queue = None
        self.__found_correct = None

    def WriteBins(self, F_HANDLE):
        """
        WriteBins
            Creates an output file describing the binning solution.
        """

        for this_elm in self.elm_order:
            F_HANDLE.write(this_elm + '\n')
            for this_bin in self.precalc_elm_bin_dict[this_elm]:
                F_HANDLE.write(str(this_bin) + '\t')
            F_HANDLE.write('\n')
        





        
