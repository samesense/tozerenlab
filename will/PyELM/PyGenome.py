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




class PyGenome():
    def __init__(self, REF_SEQ, TEST_SEQ):

        self.test_seq = TEST_SEQ
        self.ref_seq = REF_SEQ
        self.alignment = None
        self.mapping = None


    def AlignSeqs(self):
        """
        AlignSeqs
            Aligns the TEST_SEQ to the REF_SEQ using the PyAlign clustalW
        interface.  It then creats a mapping standardizes the test sequence to
        the reference provided reference sequence.

        """

        control = PyAlign.ClustalInterface(INPUT_SEQS = [self.ref_seq,
                                                         self.test_seq])

        self.alignment = control.AlignSeqs()

    
