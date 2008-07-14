from __future__ import with_statement
import PyVirus
import shelve
import numpy
from Bio import SeqIO
import itertools
import copy
import re
import os
import threading
import copy

ref_source = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
dest_dir = "C:\\local_blast\\PyELMData\\"
bkg_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'


def MakeBatch(ITERABLE, NUM_PER_BATCH):
    """
    MakeBatch
        Zips an iterable into batches and then yields the batches
    """
    internal_list = []
    for this_iter in ITERABLE:
        internal_list.append(this_iter)
        if len(internal_list) == NUM_PER_BATCH:
            yield copy.deepcopy(internal_list)
            internal_list = []
    yield copy.deepcopy(internal_list)


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
    def __init__(self, REF_SOURCE, BKG_FILE, DEST_DIR, SHELF_NAME):
        self.ref_source = REF_SOURCE
        self.bkg_source = BKG_FILE
        self.dest_dir = DEST_DIR
        self.ref_base = PyVirus.RefBase(REF_SOURCE, DEST_DIR)
        self.my_shelf_name = DEST_DIR + SHELF_NAME
        self.my_shelf = shelve.open(self.my_shelf_name)

    def BuildRefBase(self):
        """
        BuildRefBase
            Build the reference BLAST databases.
        """
        self.ref_base.BuildDatabase()

    def CreateShelf(self, MULTI_THREAD = False):
        """
        CreateShelf
            Blasts the background sequences against the RefBase and creates
        shelf instances of the resulting alignments.
        """

        BATCH_SIZE = 1
        
        with open(self.bkg_source, mode = 'r') as bkg_handle:
            bkg_iter = SeqIO.parse(bkg_handle, 'fasta')
        
            if not(MULTI_THREAD):    
                for this_seq_batch in MakeBatch(bkg_iter, BATCH_SIZE):
                    this_blast = self.ref_base.BLASTn(this_seq_batch)
                    this_seq_batch.reverse()
                    
                    for this_mapping in self.MapToRefs(this_blast,
                                                       this_seq_batch):
                        volume_name = this_mapping.ref_name
                        volume_name += this_mapping.test_viral_seq.seq_name

                        self.my_shelf[volume_name] = this_mapping

                
                
                
        
    def MapToRefs(self, BLAST_RECORDS, SEQ_RECORD):
        """
        MapToRefs
            Use a sequence of BLAST_RECORDs to create a sequence of
            MappingRecord.
        """
        #reg_exp gets the Genbank ID from the header

        for this_rec in itertools.izip(BLAST_RECORDS, iter(SEQ_RECORD)):
            checker = re.match('.*?\|(\w{1,2}\d*\.\d*).*', this_rec[0].title)
            ref_name = str(checker.groups()[0])
            this_bkgseq = PyVirus.BkgSeq(this_rec[1], None)
            this_mapping = MappingRecord(ref_name, this_bkgseq)
            this_mapping.CalculateMapping(this_rec[0])

            yield this_mapping

    
        
        
