from __future__ import with_statement
import optparse, sys, os, string
import HIVDatabase, PyVirus, AnnotUtils
from Bio import SeqIO





def Fasta2ViralSeq(INPUT_FILE):
	"""
	Reads the fasta-file provided and converts them to ViralSeq objects
	for further processing.
	
	@param INPUT_FILE:		A fasta formatted input file of nucleotides
	
	@returns				A list of ViralSeq objects
	"""
	
	output_list = []
	
	with open(INPUT_FILE) as handle:
		for this_seq in SeqIO.parse(handle, 'fasta'):
			output_list.append(PyVirus.ViralSeq(this_seq, None))
	
	return output_list
	
def ParseOptions(INPUT_STRING):
	"""
	Parses the input options
	
	"""
	
	
	
	parser = optparse.OptionParser()
	
	parser.add_option('-d', '--direc',
						dest = 'cache_dir',
						type = 'string',
						help = 'A working directory for file storage')
	parser.add_option('-o', '--out_direc',
						dest = 'out_direc',
						type = 'string',
						default = os.curdir,
						help = 'An output directory')
	parser.add_option('-n', '--num_threads',
						dest = 'num_threads',
						type = 'int',
						default = 10,
						help = 'Number of threads to process')
	
	#take input sequences
	seq_group = optparse.OptionGroup(parser, 'Seqeunce Inputs')
	
	seq_group.add_option('-i', '--input_fasta',
						dest = 'input_fasta',
						type = 'string',
						action = 'append',
						help = 'Whole Genome Sequences to visualize')
	seq_group.add_option('--ref_direc',
						dest = 'ref_direc',
						type = 'string',
						help = 'Directory of Reference Sequences')
	
	parser.add_option_group(seq_group)
	
	#options for types of outputs
	feat_group = optparse.OptionGroup(parser, 'Feature Outputs')
	feat_group.add_option('--ELM',
							dest = 'elm_flag',
							action = 'store_true',
							default = False,
							help = 'Toggle ELM Features')
	feat_group.add_option('--TF',
							dest = 'TF_flag',
							action = 'store_true',
							default = False,
							help = 'Toggle TF Features')
	feat_group.add_option('--MIRNA',
							dest = 'mirna_flag',
							action = 'store_true',
							default = False,
							help = 'Toggle miRNA Features')
	feat_group.add_option('--HOM',
							dest = 'hom_flag',
							action = 'store_true',
							default = False,
							help = 'Toggle Homology Islands')
	feat_group.add_option('--KownHom',
							dest = 'hom_file',
							type = 'string',
							help = 'A file of known Homology Islands')
	
	parser.add_option_group(feat_group)
	
	#options for defining externa dependancies
	dep_group = optparse.OptionGroup(parser, 'Dependency Locations')
	
	dep_group.add_option('--blast_dir',
						dest = 'blast_dir',
						default = string.join(['C:', 'local_blast'], 
												os.sep),
						type = 'string',
						help = 'Location of Local Blast binaries')
	dep_group.add_option('--hybrid_dir',
						dest = 'hybrid_dir',
						default = string.join(['C:', 'RNAHybrid'], 
												os.sep),
						type = 'string',
						help = 'Location of Hybrid binaries')
	dep_group.add_option('--hybrid_calib',
						dest = 'hybrid_calib',
						default = string.join([os.environ['PYTHONSCRATCH'],
												'RNAiCalibrations_KEEP.pkl'], 
												os.sep),
						type = 'string',
						help = 'Location of Hybrid binaries')
	dep_group.add_option('--elm_direc',
						dest = 'elm_direc'
						default = string.join([os.environ['MYDOCPATH'],
									'ELM_Motif_finder',
									'ELM_RAW_DOWNLOAD'], os.sep)
						type = 'string',
						help = 'Location of ELM website Files')
	
	
	parser.add_option_group(dep_group)
	
	
	
	(options, args) = parser.parse_args(INPUT_STRING)
	
	return options, args


	
class VisualController():
	def __init__(self, OPTIONS, ARGS):
		"""
		Initializes the controller class using the output created by optparse.
		"""
		
		self.options = OPTIONS
		self.args = ARGS
	
	def ProcessOptions(self):
		"""
		Processes the options provided
		"""
		
		#process input sequences
		self.test_seqs = []
		for this_file in self.options.input_fasta:
			self.test_seqs += Fasta2ViralSeq(this_file)
		
		ref_files = self.option.ref_files
		
		#process reference sequences
		
		self.mapping_base = PyVirus.HIVDatabase(self.options.cache_dir, 
												self.options.ref_direc,
												'SEQ_FILE')
		
		self.mapping_base.BuildRefBase()
		self.AddtoShelf(self.test_seqs)
		
	def run(self):
		"""
		Performs actual computation after loading options
		"""
		
		ref_base = self.mapping_base.ref_base
		
		if self.options.elm_flag:
			elm_direc = AnnotUtils.ELMParser(DIRECTORY = self.options.elm_direc)
		if self.options.mirna_flag:
			
		
		
		for this_seq_name in self.mapping_base.test_names:
			this_seq = self.mapping_base.my_test_shelf[this_seq_name]
			if self.options.elm_flag:
				this_seq.TranslateAll(ref_base)
				this_seq.FindELMs(elm_direc)
			if self.options.tf_flag:
				this_seq.FindTFSites()
			if self.options.mirna_flag:
				this_seq.HumanMiRNAsite()
				
	





if __name__ == '__main__':
	
	options, args = ParseOptions(sys.argv[1:])
	
	this_visual = VisualController(options, args)
	
	this_visual.run()
	
	#start doing stuff with input options
	
	
	
	
	
	
	
	
	
	


