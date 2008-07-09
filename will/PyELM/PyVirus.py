import shelve
import PyAlign
import itertools as IT


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
                raise NotImplemented

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
        def __init__(self, TEST_SEQ):
                self.my_sequence = TEST_SEQ
                self.known_subtype = None
                self.seq_name = None
                self.protein_starts = None
                self.protein_seqs = None














        
        


    
