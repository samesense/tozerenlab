import nose.tools
import os
import WebLogoGenerate


def testWebLogo():
	"""
	Test the weblogo creation from a list of sequences
	"""
	
	test_seqs = [	'ACTCTCATGAT',
					'ACTCTGATGAT',
					'ACTCTGATGAT',
					'ACTCTTATGAT',
					'ACTCTGATGAT',
					'ACTCTCATGAT',
					'ACTCTTATGAT',
					'ACTCTAATGAT',
					'ACTCTAATGAT',
					'ACTCTAATGAT',
					'ACTCTCATGAT']
	out_file = os.environ['MYDOCPATH'] + 'HomologyPage\\test_fig.eps'
	
	WebLogoGenerate.GenerateSeqLogo(test_seqs, out_file)
	
	os.remove(out_file)