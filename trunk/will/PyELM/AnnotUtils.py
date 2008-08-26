from __future__ import with_statement
from Bio import SeqFeature as SeqFeatClass
#a hack but the only way I can get everything to work right
from Bio.SeqFeature import SeqFeature
from reportlab.lib import colors
import pickle
import itertools
import subprocess
import re
import os
import PyMozilla
import string
import copy

class Annot():
	def __init__(self, NAME, START_POS, END_POS, HOM):
		self.name = NAME
		self.start = START_POS
		self.end = END_POS
		self.type = 'default'
		self.color = colors.black
		self.hom = HOM
		
	def __str__(self):
		return self.type + ":" + self.name + ":" + str(self.start)
	
	def GetSeqFeature(self):
		"""
		Returns a BioPython SeqFeature describing the gene represented
		"""
		seq_loc = SeqFeatClass.FeatureLocation(self.start, self.end)
		seq_feat = SeqFeature(location = seq_loc, strand = 1, 
								id = self.name)
		return seq_feat
	
	def GeneAnnot(self, GENE, REL_START, REL_STOP):
		"""
		Annotate the gene position where this ELM occurs
		"""
		self.gene = GENE
		self.rel_start = REL_START
		self.rel_stop = REL_STOP
	
	def GetSeqFeature_PROT(self):
		"""
		Returns a BioPython SeqFeature describing the ELM represented
		in reference to the protein it is on.
		"""
		seq_loc = SeqFeatClass.FeatureLocation(self.rel_start, self.rel_stop)
		seq_feat = SeqFeature(location = seq_loc, strand = 1, 
								id = self.name)
		return seq_feat
		
	def MapMeToUS(self, EQ, ANCHOR):
		"""
		Returns a new instance of self in which the start-end coordinates are 
		an average mapping to the provided coordinates
		
		EQ			Equivelent features on THIS GENOME
		ANCHOR		An ANCHORED set of features
		"""
		
		output_mapping = copy.deepcopy(self)
		
		disp = 0
		count = 0
		#for this_mapping in zip(EQ, ANCHOR):
		this_mapping = (EQ, ANCHOR)
		if (this_mapping[0] != None) & (this_mapping[1] != None):
			disp = this_mapping[0].start - this_mapping[1].start
			count = 1
			#break
	
		if count == 0:
			return
		disp_mean = disp/count
		
		output_mapping.start -= disp_mean
		output_mapping.end -= disp_mean
		
		if output_mapping.start < 0:
			return
		
		return output_mapping
	
	def CheckRange(self, TYPE, POS, FUDGE_POS = 1000):
		"""
		Returns True if this Annotation is of the proper type and starts 
		within the FUDGE region
		"""
		
		if self.type != TYPE:
			return False
		
		if abs(self.start - POS) < FUDGE_POS:
			return True
		else:
			return False

class Gene(Annot):
	def __init__(self, NAME, PRODUCT, NT_SEQ, AA_SEQ, START, END):
		self.nt_seq = NT_SEQ
		self.aa_seq = AA_SEQ
		self.start = START
		self.end = END
		self.name = NAME
		self.product = PRODUCT
		self.type = 'GENE'
		self.rel_start = None
		self.rel_stop = None

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
		self.hom = HOM

	def GetSeqFeature_PROT(self):
		raise NotImplemented
	def GeneAnnot(self, GENE, REL_START, REL_STOP):
		raise NotImplemented
		
class HomIsland(Annot):
	def __init__(self, SEQ, START_POS, END_POS, HOM):
		self.seq = SEQ
		self.start = START_POS
		self.hom = HOM
		self.name = SEQ
		self.end = END_POS
		self.type = 'HomIsland'
		self.color = colors.green
	def GetSeqFeature_PROT(self):
		raise NotImplemented
	def GeneAnnot(self, GENE, REL_START, REL_STOP):
		raise NotImplemented
	
	def __hash__(self):
		return hash(self.seq)
		
class ELM(Annot):
	def __init__(self, NAME, START_POS, END_POS, HOM):
		self.start = START_POS
		self.hom = HOM
		self.name = NAME
		self.end = END_POS
		self.type = 'ELM'
		self.gene = None
		self.rel_start = None
		self.rel_stop = None
		
		if NAME[0:3] == 'CLV':
			self.color = colors.red
		elif NAME[0:3] == 'LIG':
			self.color = colors.purple
		elif NAME[0:3] == 'MOD':
			self.color = colors.brown
		elif NAME[0:3] == 'TRG':
			self.color = colors.pink
		else:
			self.color = colors.black
	

class TFSite(Annot):
	def __init__(self, NAME, START_POS, END_POS, HOM):
		self.start = START_POS
		self.hom = HOM
		self.name = NAME
		self.end = END_POS
		self.type = 'TFSite'
		self.color = colors.silver
	
	def GetSeqFeature_PROT(self):
		raise NotImplemented
	def GeneAnnot(self, GENE, REL_START, REL_STOP):
		raise NotImplemented
		
		
		
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
	
	if len(CHROM_SEQ) > 600:
		#sequence is too long so it must be split up so recursivel call the 
		#function with smaller OVERLAPPING segments and append them together
		step_fun = lambda x:(max(0,(x)*500-100), (x+1)*500)
		off_set_fun = lambda x,y: [x[0]+y,x[1]+y,x[2]+y]
		output_list = []
		for this_win in itertools.imap(step_fun, itertools.count()):
			if this_win[1] < len(CHROM_SEQ):
				this_win_seq = CHROM_SEQ[this_win[0]:this_win[1]]
				break_val = False
			else:
				this_win_seq = CHROM_SEQ[this_win[0]:]
				break_val = True
			this_output = HybridSeq(RNA_SEQ, DELTA_THETA, this_win_seq)
			if len(this_output) > 0:
				this_output = list(itertools.imap(off_set_fun, 
												iter(this_output),
												itertools.repeat(this_win[0])))
				output_list += this_output
			if break_val:
				break
		return output_list
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
		convert_fun = lambda x: [int(x[2]), float(x[0]), float(x[1])]
		final_output = map(convert_fun, final_output)
		
	return final_output
	
def ELMParser(DIRECTORY = None):
	"""
		Parser
			Parses a directory of HTML files to extract the ELM names and
			Regular Expressions and instanciates it into SELF.

			Parser(DIRECTORY)

			DIRECTORY		The directory containing the HTML files downloaded
							from the ELM database.  If left empty then
			C:\Documents and Settings\Will\My Documents\PyELM\ELM_RAW_DOWNLOAD\

	"""
	if DIRECTORY == None:
		DIRECTORY = os.environ['MYDOCPATH'] + "ELM_Motif_finder\\ELM_RAW_DOWNLOAD\\"

	file_list = os.listdir(DIRECTORY)
	all_reg_exps = {}

	for i in xrange(1, len(file_list) - 1):
		this_file = open(DIRECTORY + file_list[i], 'r')
		elm_name = file_list[i].rpartition('.')
		current_line = ''
		while current_line.find('Pattern:') == -1:
			current_line = this_file.readline()
		current_line = this_file.readline()
		
		reg_exp = current_line[current_line.find('>') +
												 1:current_line.find('/') - 2]
		comp_reg_exp = re.compile(reg_exp)
		all_reg_exps[elm_name[0]] = (reg_exp, comp_reg_exp)
		this_file.close()
	
	return all_reg_exps
	
	
def TFChecker(INPUT_SEQ):
	"""
	Uses a MATCH login to find the transcription factor binding sites in 
	INPUT_SEQ.
	
	Returns a list of lists:
		[[pos1, strand1, core_sim1, total_sim1, seq1, name1, mat1],
		 [[pos2, strand2, core_sim2, total_sim2, seq2, name2, mat2]]
	
	"""
	moz_emu = PyMozilla.MozillaEmulator(cacher=None)
	info_getter = re.compile('<A HREF.*?</a>(.*?)<a href.*?</a>',re.S)
	mat_getter = re.compile('<A HREF.*?>(.*?)</a>')
	name_getter = re.compile('<a href.*?>(.*?)</a>')
	sep = re.compile(' *')
	base_url = 'http://www.gene-regulation.com'
	user_name = 'walldo2'
	password = 'bfg20k'
	login_post = 'user=' + user_name
	login_post += '&password=' + password
	login_post += '&request_uri=%2Findex.html&Log+In=Log+In'
	
	#login and get the cookie
	moz_emu.download(base_url + '/login', login_post)
	
	match_base = '/cgi-bin/pub/programs/match/bin/match.cgi'
	post_data = 'Status=First'
	post_data += '&searchName=default'
	post_data += '&usr_seq=default.seq'
	post_data += '&seqStat=DEL'
	post_data += '&sequenceName=default'
	post_data += '&&theSequence=' + INPUT_SEQ
	post_data += '&SearchMode=ALL'
	post_data += '&group=vertebrates'
	post_data += '&quality=high'
	post_data += '&allMode=FP'
	post_data += '&cut-off=0.7'
	post_data += '&core=0.75'
	post_data += '&ourProfile=muscle_specific.prf'
	
	output = moz_emu.download(base_url + match_base,post_data)
	
	val = map(string.strip,info_getter.findall(output))
	info = map(sep.split,val)
	tf_names = name_getter.findall(output)
	mat_names = mat_getter.findall(output)
	
	if len(val) == 0:
		return []
	
	converter_fun = lambda x:[int(x[0]), x[1], float(x[2]), float(x[3]), x[4]]
	
	final_output = []
	for i in range(len(info)):
		this_out = converter_fun(info[i])
		this_out.append(tf_names[i])
		this_out.append(mat_names[i])
		final_output.append(this_out)
	return final_output