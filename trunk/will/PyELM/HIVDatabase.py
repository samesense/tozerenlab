from __future__ import with_statement
import PyVirus
import numpy
from Bio import SeqIO
import itertools
import copy
import re
import os
import threading
import copy
import logging
import cPickle as pickle
import types

ref_source = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
dest_dir = "C:\\local_blast\\PyELMData\\"
bkg_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

class MappingRecord():
    def __init__(self, REF_NAME, TEST_VIRAL):
        self.ref_name = REF_NAME

        self.test_viral_seq = TEST_VIRAL

        self.mapping = None
        self.is_match = None

    def CalculateMapping(self, BLAST_ALIGNMENT):
        self.mapping = numpy.zeros((1,len(self.test_viral_seq.my_sequence)))
        self.is_match = numpy.zeros((1,len(self.test_viral_seq.my_sequence)))
        filter_fun = lambda x: x[0] == '|'
        unzip = lambda x: x[1]
        for this_align in BLAST_ALIGNMENT.hsps:
            if this_align.query_start > this_align.query_end:
                continue
            if this_align.sbjct_start > this_align.query_end:
                continue
            
            query_inds = range(this_align.query_start, this_align.query_end)
            subjct_inds = range(this_align.sbjct_start, this_align.sbjct_end)
            self.mapping[0,query_inds] = numpy.array(subjct_inds)
            
            matched_inds = filter(filter_fun, zip(this_align.match, query_inds))
            matched_inds = map(unzip, matched_inds)
            self.is_match[0,matched_inds] = 1
            

    def __hash__(self):
        return hash(self.ref_name + self.test_viral_seq.seq_name)



class MappingBase():
    def __init__(self, REF_SOURCE, DEST_DIR, SHELF_NAME):
        self.ref_source = REF_SOURCE
        self.dest_dir = DEST_DIR
        self.ref_base = PyVirus.RefBase(REF_SOURCE, DEST_DIR)
        self.my_shelf_name = DEST_DIR + SHELF_NAME
        try:
            with open(self.my_shelf_name) as prev_handle:
                self.my_shelf = pickle.load(prev_handle)
        except IOError:
            self.my_shelf = {}
            
           
        

    def BuildRefBase(self):
        """
        BuildRefBase
            Build the reference BLAST databases.
        """
        self.ref_base.BuildDatabase()

    def AddtoShelf(self, INPUT_SEQRECORD_ITER, MULTI_THREAD = False):
        """
        AddtoShelf
            Blasts the background sequences against the RefBase and creates
        shelf instances of the resulting alignments.  Input must be a LIST,
        ITER or GENERATOR
        """
        
           
        for this_seq in INPUT_SEQRECORD_ITER:
            this_blast = self.ref_base.BLASTn(this_seq)

            for this_mapping in self.MapToRefs(this_blast,
                                           this_seq):
                volume_name = this_mapping.ref_name
                volume_name += this_mapping.test_viral_seq.seq_name
                print volume_name

                self.my_shelf[volume_name] = this_mapping

    def SaveShelf(self):
        """
        Saves the precalculated dictionary
        """
        with open(self.my_shelf_name, mode = 'w') as prev_handle:
            pickle.dump(self.my_shelf, prev_handle)

        
             

    def MapToRefs(self, BLAST_RECORDS, SEQ_RECORD):
        """
            MapToRefs
            Use a sequence of BLAST_RECORDs to create a sequence of
        MappingRecord.
        """
        #reg_exp gets the Genbank ID from the header

        checker = re.compile('.*?\|(\w{1,2}\d*\.\d*).*')


        for this_align in BLAST_RECORDS.alignments:
            temp_name = checker.match(this_align.title)
            ref_name = str(temp_name.groups()[0])
            this_bkgseq = PyVirus.BkgSeq(SEQ_RECORD, None)
            this_mapping = MappingRecord(ref_name, this_bkgseq)
            this_mapping.CalculateMapping(this_align)

            yield this_mapping



