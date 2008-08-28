from __future__ import with_statement
import optparse, sys, os, string, pickle, logging
import HIVDatabase, PyVirus, AnnotUtils
from Bio import SeqIO





def Fasta2ViralSeq(INPUT_FILE, WANTED_SUBS):
	"""
	Reads the fasta-file provided and converts them to ViralSeq objects
	for further processing.
	
	@param INPUT_FILE:		A fasta formatted input file of nucleotides
	
	@returns				A list of ViralSeq objects
	"""
	
	output_list = []
	
	with open(INPUT_FILE) as handle:
		for this_seq in SeqIO.parse(handle, 'fasta'):
			this_virus = PyVirus.ViralSeq(this_seq, None)
			if WANTED_SUBS == None:
				output_list.append(this_virus)
			elif this_virus.tested_subtype in WANTED_SUBS:
				output_list.append(this_virus)
	
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
	parser.add_option('-v', '--verbose',
						dest = 'verbose',
						action = 'store_true',
						default = False,
						help = 'Be Verbose')
	parser.add_option('--log_file',
						dest = 'log_name',
						type = 'string',
						help = 'Log File Name')
	
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
						default = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\',
						help = 'Directory of Reference Sequences')
	seq_group.add_option('-s','--subtype',
						dest = 'wanted_subs',
						action = 'append',
			help = 'The desired subtypes in the output: can be used multiple times')
	
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
						default = os.environ['PYTHONSCRATCH'] + 'RNAiCalibrations_KEEP.pkl',
						type = 'string',
						help = 'Location of Hybrid binaries')
	dep_group.add_option('--elm_direc',
						dest = 'elm_direc',
						default = os.environ['MYDOCPATH'] + string.join(['ELM_Motif_finder',
									'ELM_RAW_DOWNLOAD'], os.sep) + os.sep,
						type = 'string',
						help = 'Location of ELM website Files')
	dep_group.add_option('--scratch_dir',
						dest = 'scratch_dir',
						default = string.join(['C:', 'new_scratch'], 
												os.sep),
						type = 'string',
						help = 'A scratch directory needed for data storage.')
	
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
		
		self.ProcessOptions()
	
	def ProcessOptions(self):
		"""
		Processes the options provided
		"""
		
		logging.info('Reading sequences')
		#process input sequences
		self.test_seqs = []
		for this_file in self.options.input_fasta:
			self.test_seqs += Fasta2ViralSeq(this_file, self.options.wanted_subs)
		logging.info('%(num)d sequences read' % {'num': len(self.test_seqs)})
		
		#process reference sequences
		
		logging.info('Creating MappingBase')
		self.mapping_base = HIVDatabase.MappingBase(self.options.ref_direc,
												self.options.scratch_dir, 
												'SEQ_FILE')
		
		logging.info('Building Fasta Database')
		self.mapping_base.BuildRefBase()
		
		logging.info('Adding sequences to the Shelf')
		self.mapping_base.AddtoShelf(self.test_seqs)
		
	def run(self):
		"""
		Performs actual computation after loading options
		"""
		
		ref_base = self.mapping_base.ref_base
		calib_dict = None
		elm_dict = None
		
		if self.options.elm_flag:
			elm_dict = AnnotUtils.ELMParser(DIRECTORY = self.options.elm_direc)
		if self.options.mirna_flag:
			with open(self.options.hybrid_calib) as handle:
				calib_dict = pickle.load(handle)
		for this_calib in calib_dict.keys()[10:]:
			junk = calib_dict.pop(this_calib)
		
		logging.debug('Beginning Annotation')
		self.mapping_base.AnnotateBase(calib_dict, elm_dict, True, 
									ref_base.ref_seqs[0].seq_name)
		
		(gene_fig, prot_fig) = self.mapping_base.MakeMultiDiagram(ref_base.ref_seqs[0].seq_name,
															self.options.wanted_subs)
		gene_fig.draw(format = 'linear', fragments = 1)
		gene_fig.write(self.options.out_direc + 'test_out.pdf', 'PDF')





if __name__ == '__main__':
	
	options, args = ParseOptions(sys.argv[1:])
	
	file_fmt = '%(levelname)s \t %(threadName)s \t %(funcName)s \t %(asctime)s \t %(message)s'
	console_fmt = '%(asctime)-7s %(message)s'
	
	
	
	
	if options.verbose:
		logging.basicConfig(level = logging.DEBUG, fmt = console_fmt)
	else:
		logging.basicConfig(level = logging.WARNING, fmt = console_fmt)
	
	if options.log_name != None:
		file_name = options.out_direc + options.log_name
		file_log = logging.StreamHandler(open(file_name, mode = 'w'))
		formatter = logging.Formatter(file_fmt)
		file_log.setLevel(logging.DEBUG)
		file_log.setFormatter(formatter)
		logging.getLogger('').addHandler(file_log)
	
	logging.info('Finished Input Parsing')
	
	
	this_visual = VisualController(options, args)
	
	logging.info('Finished Pre-Processing')
	
	
	this_visual.run()
	
	logging.info('Finished Running')
	
	#start doing stuff with input options
	
	
	
	
	
	
	
	
	
	


