from __future__ import with_statement
import nose.tools, os
import VisualGenome




def testParsingOptions():
	"""
	Ensures that the OptionParsing stays consisent
	"""
	
	this_dir = os.environ['MYDOCPATH'] + 'PyELM\\'
	
	poss_opt = [['-d', '--direc', 'cache_dir', 'C:\\pycratch', None],
				['-i', '--input_fasta', 'input_fasta', '50_seqs.fasta', None],
				['-r', '--ref_seqs', 'ref_files', 'Globalaf033819.gb', None],
				['-b', '--ref_direc', 'ref_direc', this_dir, None],
				[None, '--blast_dir', 'blast_dir', 'C:\\pycratch', None]]
	
	for this_opt in poss_opt:
		yield CheckOptions, this_opt
	
def CheckOptions(INPUT):
	"""
	Checks an option provided in the form:
	[short_form, long_form, dest, input_value, result]
		
		If result is None then the desired value is assumed to equal the 
		input_value
		
	"""
	
	if INPUT[4] == None:
		correct_val = INPUT[3]
	else:
		correct_val = INPUT[4]
	
	if INPUT[0] != None:
		short_str = 'python VisualGenome.py ' + INPUT[0] + ' ' + INPUT[3]
		print short_str
		(options, args) = VisualGenome.ParseOptions(short_str)
		try:
			val = options.__getattribute__[INPUT[2]]
		except AttributeError:
			nose.tools.assert_true(False, 'Missing Attribute: ' + INPUT[2])
		
		nose.tools.assert_true(val == correct_val, 
								'The attribute did not save properly')
	
	if INPUT[1] != None:
		short_str = 'python VisualGenome.py ' + INPUT[0] + '=' + INPUT[3]
		(options, args) = VisualGenome.ParseOptions(short_str)
		try:
			val = options.__getattribute__[INPUT[2]]
		except AttributeError:
			nose.tools.assert_true(False, 'Missing Attribute: ' + INPUT[2])
		
		nose.tools.assert_true(val == correct_val, 
								'The attribute did not save properly')