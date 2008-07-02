import os
import numpy
import re
from pyQueueUtils import *
from copy import *
import threading
import time
import bisect
import Queue

class PyELM:

    def __init__(self):

        self.elm_dict = {}
        self.elm_order = []
        self.sequences = []
        self.elm_re_dict = {}
        self.max_size = 0

        self.__found_correct = False

        self.precalc_hist_array = []        
        self.precalc_elm_match_seqs = []
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
        self.precalc_elm_match_seqs = []
        self.precalc_hist_array = []
        
    def LoadSeqs(self, SEQUENCES):
        """
        LoadSeqs
            Loads a set of sequences into this instance for further analysis.

            LoadSeqs(SEQUENCES)

            SEQUENCES       A LIST of characters representing Amino Acids


        """

        self.sequences = SEQUENCES
        self.max_size = 0
        for i in xrange(len(self.sequences)):
            if len(self.sequences[i]) > self.max_size :
                self.max_size = len(self.sequences[i])


    

    def ELMMatchSeqs(self):
        """
        ELMMatchSeqs
            Searches for the regular expressions of ELM_DICT within the
            SEQUENCES provided.

        ELM_MATCH_VEC = ELMMatchSeqs()


        ELM_MATCH_VEC   A LIST corresponding to each BioSeq object provided.
                        Each element in the list is a DICT where each element
                        indicates a matched ELM location.

        """

        if len(self.elm_dict) == 0:
            print 'ELMParser must be run first!'
            return
        if len(self.sequences) == 0:
            print 'LoadSeqs must be run first!'
            return


        if len(self.precalc_elm_match_seqs) == 0:
            #compile all of the re's
            if self.elm_re_dict != None:
                elm_re_dict = {}
                for this_elm in self.elm_order:
                    elm_re_dict[this_elm] = re.compile(self.elm_dict[this_elm],
                                                       re.I)
                self.elm_re_dict = elm_re_dict
            
            seq_elm_dict = []
            for this_seq in self.sequences:
                this_elm_seq = {}
                for this_elm in self.elm_order:
                    this_index = []
                    #each search only finds one match
                    #so repeat the search until no more are found
                    spot = self.elm_re_dict[this_elm].search(this_seq)
                    while spot != None:
                        this_index.append(spot.start())
                        spot = self.elm_re_dict[this_elm].search(this_seq,
                                                                 spot.start()
                                                                 + 1)
                    this_elm_seq[this_elm] = this_index
                seq_elm_dict.append(this_elm_seq)


            self.precalc_elm_match_seqs = seq_elm_dict

        
        return self.precalc_elm_match_seqs
        
    def ELMHistGen(self):
        """
        ELMHistGen
            Generates a count-matrix of ELM matches at specific locations in
            the sequences.

        ELM_COUNT_MAT = ELMHistGen()

        ELM_COUNT_MAT   An [maxSeqLength x len(ELM_DICT)] each Row represent
                        an ELM and each Column is a seqLocation.

        """
        if len(self.elm_dict) == 0:
            print 'ELMParser must be run first!'
            return
        if len(self.sequences) == 0:
            print 'LoadSeqs must be run first!'
            return

        if len(self.precalc_elm_match_seqs) == 0:
            self.ELMMatchSeqs()

        if len(self.precalc_hist_array) == 0:

            hist_array = numpy.zeros((len(self.elm_order), self.max_size))

            counter = 0
            for this_elm in self.elm_order:
                for i in xrange(len(self.precalc_elm_match_seqs)):
                    hist_array[counter,
                               self.precalc_elm_match_seqs[i][this_elm]] += 1
                counter += 1

            self.precalc_hist_array = hist_array
        
        return self.precalc_hist_array

    def SimpleELMCall(self, SEQ_IND, WANTED_ELM):
        """
        SimpleELMCall
            Returns a simple Presence/Absence Call for the ELM.

        PACall = SimpleELMCall(SEQ_IND,WANTED_ELM)

        SEQ_IND     The index of the desired sequence
        WANTED_ELM  The desired ELM

        PACall      Returns 1 if the ELM is present in the sequence and
                    0 otherwise.

        """

        if SEQ_IND > len(self.sequences) | SEQ_IND < 0:
            print 'Arguement SEQ_IND must be between 0 and len(self.sequences)'
            raise IndexError
        
        if len(self.precalc_elm_match_seqs) == 0:
            self.ELMMatchSeqs()

        if len(self.precalc_elm_match_seqs[SEQ_IND][WANTED_ELM]) >= 1:
            return 1
        else:
            return 0



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
            if sum(self.ELMPosGen(SEQ_IND, allBins[i], allBins[i+1],
                                  WANTED_ELM)) > 0:
                thisCall.append(1)
            else:
                thisCall.append(0)
        return thisCall


        


    def ELMPosGen(self, SEQ_IND, POS_START, POS_STOP, WANTED_ELM):
        """
        ELMPosGen
            Generates a [len(ELM_MATCH_VEC) MAX_SIZE] matrix in which each
            position is a 1 if the ELM is present at that location and 0
            otherwise.

            ELM_POS_DICT = ELMPosGen(POS_START,POS_STOP,WANTED_ELM)

        ELM_POS_DICT    A DICT which is 'keyed' by ELM name and each value
                        is a matrix.


            ELM_POS_ARRAY = ELMPosGen(..., WANTED_ELM)
        Will return only the desired array.

        """

        #check to see if calculation can be performed
        if len(self.precalc_elm_match_seqs) == 0:
            self.ELMMatchSeqs()

        if POS_START < 0 | POS_STOP <= POS_START | POS_START > self.max_size | POS_STOP > self.max_size:
            print 'Arguements to POS_START >0 & POS_STOP > POS_START'
            raise IndexError

        if SEQ_IND > len(self.sequences) | SEQ_IND < 0:
            print 'Arguement SEQ_IND must be between 0 and len(self.sequences)'
            raise IndexError

        if not(WANTED_ELM in self.elm_order):
            print 'WANTED_ELM not found in elm_order'
            raise KeyError

        boolArray = numpy.zeros((1, POS_STOP - POS_START), 'bool')
        for i in xrange(POS_START, POS_STOP):
            boolArray[0, i - POS_START] = i in self.precalc_elm_match_seqs[SEQ_IND][WANTED_ELM]
            
        return boolArray

    def GetELMIterator(self):
        """
        GetELMIterator
            Returns an iterator which will iterate over the items in
            the elm_dict.
        """
        return self.ELMIterator(self)
        

    class ELMIterator():
        """
        ELMIterator
            An interator which passes over multiple ELMS
        """
        def __init__(self, INSTANCE):
            self.__list = INSTANCE.elm_order
            self.__this_spot = -1

            
        def __iter__(self):
            return self

        def next(self):
            """
            next
                Procededs to the next ELM
            """
            self.__this_spot += 1
            if self.__this_spot < len(self.__list):
                return self.__list[self.__this_spot]
            else:
                raise StopIteration
    
    def GetSeqIterator(self, WANTED_ELM = None, WANTED_RANGE = None):
        """
        GetSeqIterator
            Will allow for easy iteration and retrival of data in a
            Seq-Order manner.

        GetSeqIterator('ELMPosDict',WANTED_ELM,WANTED_RANGE)
            WANTED_RANGE        A tuple indicating the range of positions
                                desired.  If left empty then the entire
                                sequence is provided.
            
            This creates an iterator which will return the ELMPos Array data
            in Sequence order. If WANTED_ELM is left empty then Data is
            returned as an Array in which each row represents one ELM.
        """
        if (WANTED_ELM == None) | (WANTED_ELM in self.elm_order):
            if WANTED_RANGE == None:
                WANTED_RANGE = (0, self.max_size - 1)
            elif len(WANTED_RANGE) != 2:
                print 'WANTED_RANGE must be len == 2 or None'
                raise IndexError
            elif WANTED_RANGE[0] < 0 | WANTED_RANGE[1] <= WANTED_RANGE[0]:
                print 'WANTED_RANGE[0] > 0 and WANTED_RANGE[1] > WANTED_RANGE[0]'
                raise IndexError
            elif WANTED_RANGE[0] > self.max_size | WANTED_RANGE[1] > self.max_size:
                print 'WANTED_RANGE[0] < max_size and WANTED_RANGE[1] < max_size'
                raise IndexError
               
            return self.SeqIterator(self, WANTED_ELM, WANTED_RANGE)
        else:
            print 'Provided an unknown ELM key'
            raise KeyError
            

    class SeqIterator:
        """
        SeqIterator
            An iterator which passes over all sequences in the database
        """
        def __init__(self, INSTANCE, WANTED_ELM, WANTED_RANGE):
            self.__wanted_elm = WANTED_ELM
            self.__wanted_range = WANTED_RANGE
            self.__list_size = len(INSTANCE.sequences)
            self.__this_spot = -1
            self.__max_size = INSTANCE.max_size
            self.__elm_order = INSTANCE.elm_order
            self.elm_order = INSTANCE.elm_order
            self.precalc_elm_match_seqs = INSTANCE.precalc_elm_match_seqs
            self.ELMPosGen = INSTANCE.ELMPosGen

        def __iter__(self):
            return self

        def next(self):
            """
            next
                Procedes onto the next sequence in the database.
            """
            self.__this_spot += 1
            if self.__this_spot < self.__list_size:
                if self.__wanted_elm != None:
                    return self.ELMPosGen(self.__this_spot,
                                          self.__wanted_range[0],
                                          self.__wanted_range[1],
                                          self.__wanted_elm)
                else:
                    out_array = numpy.zeros((len(self.__elm_order),
                                             self.__wanted_range[1] -
                                             self.__wanted_range[0]))
                    
                    for i in xrange(len(self.__elm_order)):
                        out_array[i,:] = self.ELMPosGen(self.__this_spot,
                                                        self.__wanted_range[0],
                                                        self.__wanted_range[1],
                                                        self.__elm_order[i])
                    return out_array
            else:
                raise StopIteration

    def __QueueWorker(self, WANTED_ELM):
        """
        A helper function to manage the Priority Queue and Threads
        """
        while True:
            self.__global_counter += 1 
            try:
                print 'Worker qsize: ', self.__working_queue.qsize()
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
        for this_seq in self.GetSeqIterator(WANTED_ELM,
                                            (BIN_EDGES[0], BIN_EDGES[1])):
            current_bin = sum(this_seq)
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
        num_thread = 10
        if len(self.precalc_elm_match_seqs) == 0:
            self.ELMMatchSeqs()

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
        





        
