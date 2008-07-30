from __future__ import with_statement
from Bio import SeqFeature as SeqFeatClass
#a hack but the only way I can get everything to work right
from Bio.SeqFeature import SeqFeature
from reportlab.lib import colors

class Annot():
	def __init__(self, NAME, START_POS, END_POS, HOM):
		self.name = NAME
		self.start = START_POS
		self.end = END_POS
		self.type = 'default'
		self.color = colors.black
		
	def GetSeqFeature(self):
		"""
		Returns a BioPython SeqFeature describing the gene represented
		"""
		seq_loc = SeqFeatClass.FeatureLocation(self.start, self.end)
		seq_feat = SeqFeature(location = seq_loc, strand = 1, 
								id = self.name)
		return seq_feat

class Gene(Annot):
	def __init__(self, NAME, PRODUCT, NT_SEQ, AA_SEQ, START, END):
		self.nt_seq = NT_SEQ
		self.aa_seq = AA_SEQ
		self.start = START
		self.end = END
		self.name = NAME
		self.product = PRODUCT
		self.type = 'GENE'

	def __str__(self):
		temp_str = '(Gene Class '
		temp_str += 'Name: ' + self.name
		temp_str += ', Start: ' + str(self.start) + ')'
		return temp_str
		
	def __hash__(self):
		return hash(self.aa_seq)
		
class HumanMiRNA(Annot):
	def __init__(self, NAME, START_POS, END_POS, HOM):
		self.name = NAME
		self.start = START_POS
		self.end = END_POS
		self.type = 'HumanMiRNA'
		self.color = colors.blue
		
class HomIsland(Annot):
	def __init__(self, SEQ, START_POS, END_POS, HOM):
		self.seq = SEQ
		self.start = START_POS
		self.hom = HOM
		self.name = SEQ
		self.end = START_POS + len(SEQ)
		self.type = 'HomIsland'
		self.color = colors.green
		
	def __hash__(self):
		return hash(self.seq)
		
def CalibrateRNAi(INPUT_FILE, OUTPUT_FILE):
	mi_rna_dict = {}
	splitter = re.compile('\t')
	with open(INPUT_FILE) as file_handle:
		#get rid of comment line
		this_line = file_handle.next()
		for this_line in file_handle:
			data = splitter.split(this_line)
			mi_rna_dict[data[3]] = string.upper(string.replace(data[5], '-', ''))
	
	RNAHybridPath = 'C:\\RNAHybrid\\'
	CODING_FILE = 'C:\\RNAHybrid\\coding_HIV.txt'
	NUM_SAMPLES = 5000
	
	calib_dict = {}
	
	mi_keys = mi_rna_dict.keys()
	
	
	for this_mirna in mi_keys:
		
		this_seq = mi_rna_dict[this_mirna]
		command = RNAHybridPath + 'RNAcalibrate'
		
		command = RNAHybridPath + 'RNAcalibrate'

		command += ' -d ' + CODING_FILE
		command += ' -k ' + str(NUM_SAMPLES)

		command += ' ' + this_seq

		sys_call = subprocess.Popen(command, shell = True,
										stdout = subprocess.PIPE)

		sys_call.wait()
		output = sys_call.communicate()[0]
		outputList = re.split(' |\n', output)
		
		calib_dict[this_mirna] = (this_seq, outputList[2], outputList[3])
	
	with open(OUTPUT_FILE, mode = 'w') as handle:
		pickle.dump(calib_dict, handle)

def HybridSeq(RNA_SEQ, DELTA_THETA, CHROM_SEQ,
              RNA_HYBRID_PATH = 'C:\\RNAHybrid\\'):
	"""
	HybridSeq
		Makes a call to the RNAhybrid.exe to determine the likilyhood of the
		miRNA sequence binding to the Chomrosome Sequence.  It uses the data
		from Calibrate to determine the statistical significance.

			HybridSeq(RNA_SEQ,DELTA_THETA,CHROM_SEQ, 
				RNA_HYBRID_PATH = 'C:\\RNAHybrid\\'):

			RNA_SEQ         The sequence of the short miRNA.

			DELTA_THETA      The output of Calibrate ... the shape of the extreme
							value distribution.

	        RNA_HYBRID_PATH   The path to the folder that contains RNAhybrid


	    RETURNS

	        the list of tuples:
	            (position, energy, p-value)


	"""
	
	if len(CHROM_SEQ) > 700:
		#sequence is too long so it must be split up so recursivel call the 
		#function with smaller OVERLAPPING segments and append them together
		step_fun = lambda x:(max(0,(x)*500-100), (x+1)*500)
		off_set_fun = lambda x,y: [x[0]+y,x[1]+y,x[2]+y]
		output_list = []
		for this_win in itertools.imap(step_fun, itertools.count()):
			try:
				this_win_seq = CHROM_SEQ[this_win[0]:this_win[1]]
				break_val = False
			except:
				this_win_seq = CHROM_SEQ[this_win[0]]
				break_val = True
			
			this_output = HybridSeq(RNA_SEQ, DELTA_THETA, this_win_seq)
			if len(this_output) > 0:
				this_output = list(itertools.imap(off_set_fun, 
												iter(this_output),
												itertools.repeat(this_win[0])))
				output_list += this_output
				
	command = RNA_HYBRID_PATH + 'RNAhybrid'

	command += ' -c'
	command += ' -d ' + DELTA_THETA[0] + ',' + DELTA_THETA[1]

	command += ' -p 0.1' 
	command += ' ' + CHROM_SEQ
	command += ' ' + RNA_SEQ

	sys_call = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE)
	sys_call.wait()

	output = sys_call.communicate()[0]

	output_reg_exp = '(.*?):.*?command_line:[\-\.\d]*:([\.\-\d]*):'
	output_reg_exp += '([\.\-\d]*):([\.\-\d]*).*'

	final_output = []
	if len(output) > 1:
		all_lines = re.split('\n', output)

		for this_line in all_lines[0:-1]:
			final_output.append(re.match(output_reg_exp,
											this_line).groups()[1:])
		print final_output
		convert_fun = lambda x: [int(x[2]), float(x[0]), float(x[1])]
		final_output = map(convert_fun, final_output)
		
	return final_output