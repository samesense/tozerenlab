from __future__ import with_statement
import shelve
import PyAlign
import itertools as IT
import HIVGenoTypeChecker as GenoTyping
from Bio import SeqIO


def MakeMapping(ALIGNMENT):
	REF_ITER = iter(ALIGNMENT[1])
	TEST_ITER = iter(ALIGNMENT[0])
	filterFun = lambda x: x[0] != '-'
	mapping_val = []
	num_val = []
	ref_num_inter = IT.chain(IT.ifilter(filterFun, IT.izip(REF_ITER, IT.count())), IT.repeat(('-',len(ALIGNMENT[1]))))
	test_iter = IT.ifilter(filterFun, IT.izip(TEST_ITER, ref_num_inter))
	for this_iter in test_iter:
		mapping_val.append(this_iter[0])
		num_val.append(this_iter[1][1])
	print ALIGNMENT[1] + str(len(ALIGNMENT[1]))
	print ALIGNMENT[0]
	print str(mapping_val) + str(len(mapping_val))
	print str(num_val) + str(len(num_val))

class Gene():
	def __init__(self, NT_SEQ, AA_SEQ, START, END):
		self.nt_seq = NT_SEQ
		self.aa_seq = AA_SEQ
		self.start = START
		self.end = END


class ViralSeq():
        def __init__(self, TEST_SEQ):
                self.my_sequence = TEST_SEQ
                self.tested_subtype = None
                self.seq_name = None

        def DetSubtype(self):
                """
                DetSubtype
                        Determines the subtype of the sequence using the
                HIVSubtype module.
                """

                self.tested_subtype = GenoTyping.GetSimple(self.my_sequence)

class BkgSeq(ViralSeq):
        def __init__(self, TEST_SEQ):
                self.my_sequence = TEST_SEQ
                self.tested_subtype = None
                self.seq_name = None


class PatSeq(ViralSeq):
        def __init__(self, TEST_SEQ):
                self.my_sequence = TEST_SEQ
                self.pat_data = None
                self.tested_subtype = None
                self.seq_name = None


class RefSeq(ViralSeq):
        def __init__(self, FILENAME):
                self.my_sequence = None
                self.known_subtype = None
                self.seq_name = None
                self.annotation = None

                self.ParseGenbank(FILENAME)

                

        def ParseGenbank(self, FILENAME):
                with open(FILENAME, mode='r') as handle:
                        this_record = SeqIO.parse(handle, 'genbank').next()

                self.my_sequence = this_record.seq.tostring()
                

                self.annotation = {}
                for feat in this_record.features:
                        if feat.qualifiers.has_key('translation'):
                                start_pos = feat.location._start.position
                                end_pos = feat.location._end.position
                                trans_data = feat.qualifiers['translation'][0]
                                nt_seq = self.my_sequence[start_pos:end_pos]
                                this_gene = Gene(nt_seq, trans_data,
                                                 start_pos, end_pos)

                                self.annotation[feat.qualifiers['gene'][0]] = this_gene
	














        
        


    
