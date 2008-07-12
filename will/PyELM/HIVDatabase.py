from __future__ import with_statement
import PyVirus
import shelve
import numpy
from Bio import SeqIO
import itertools
import copy
import re
import os

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

    def CreateShelf(self, MAX_SEQS = None, MULTI_THREAD = False):
        """
        CreateShelf
            Blasts the background sequences against the RefBase and creates
        shelf instances of the resulting alignments.
        """
        with open(self.bkg_source, mode = 'r') as bkg_handle:
            bkg_iter = SeqIO.parse(bkg_handle, 'fasta')
        
            if not(MULTI_THREAD):    
                for this_seq in bkg_iter:
                    this_blast = self.ref_base.BLASTn(this_seq)
                    maprecord_iter = itertools.imap(self.MapToRefs,
                                                    iter(this_blast.alignments),
                                                    itertools.repeat(this_seq))
                    
                    for this_mapping in itertools.izip(maprecord_iter,
                                                       itertools.repeat(None,
                                                       MAX_SEQS)):
                        volume_name = this_mapping[0].ref_name
                        volume_name += this_mapping[0].test_viral_seq.seq_name
                        print volume_name

                        self.my_shelf[volume_name] = this_mapping[0]
                
                
        
    def MapToRefs(self, BLAST_RECORD, SEQ_RECORD):
        """
        MapToRefs
            Uses a BLAST_RECORD to create a MappingRecord.
        """
        #reg_exp gets the Genbank ID from the header
        checker = re.match('.*?\|(\w{1,2}\d*\.\d*).*',BLAST_RECORD.title)
        ref_name = str(checker.groups()[0])
        this_record = PyVirus.BkgSeq(SEQ_RECORD.seq.tostring(), SEQ_RECORD.id)
        this_mapping = MappingRecord(ref_name, this_record)
        this_mapping.CalculateMapping(BLAST_RECORD)

        return this_mapping

    
        
        
