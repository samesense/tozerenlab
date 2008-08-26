import nose.tools
import PyMatch
import time

def testGetSequence():
	"""
	Test the ability to gather sequence
	"""
	
	#a few of hand-tested genome positions
	test_data = [	('1',500,520,'GTCTGACCTGAGGAGAACTGT'),
					('2',500,520,'CCCGACCCCGACCCCGACCCA'),
					('3',50000,50020,'TCTTCTTTTATGAAAAAGGAT'),
					('4',50000,50020,'AGAGCCCTGCAATTTGAAGAT'),
					('5',100000,100020,'AATGTTCACCAGTATATTTTA'),
					('X',100000,100020,'TAGGTCTCATTGAGGACAGAT'),
					('Y',100000,100020,'TAGGTCTCATTGAGGACAGAT')]
					
	for this_check in test_data:
		yield CheckGetSequence, this_check
	
def CheckGetSequence(INPUT_DATA):
	"""
	Checks the data yielded by testGetSequence
	"""
	time.sleep(3)
	seq_data = PyMatch.GetSequence(INPUT_DATA[0], INPUT_DATA[1], INPUT_DATA[2])
	
	bad_str = 'CHR: %(chr)s pos: [%(start)s %(end)s]' % {'chr':INPUT_DATA[0], 
													'start':str(INPUT_DATA[1]),
													'end':str(INPUT_DATA[2])}
	bad_str += '%(known)s Retrieved sequence: %(ret)s' % {'known':INPUT_DATA[3],
													'ret':seq_data}
	
	nose.tools.assert_true(seq_data == INPUT_DATA[3], bad_str)