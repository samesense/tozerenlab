from __future__ import with_statement
import optparse, sys, os, string, pickle, logging, ConfigParser, types
import HIVDatabase, PyVirus, AnnotUtils, PatUtils
import pylab

from Bio import SeqIO
from collections import defaultdict




def Fasta2ViralSeq(INPUT_FILE, WANTED_SUBS, PICT_FILE):
	"""
	Reads the fasta-file provided and converts them to ViralSeq objects
	for further processing.
	
	@param INPUT_FILE:		A fasta formatted input file of nucleotides
	@param WANTED_SUBS:		A list of subtypes that are to be processed.  If 
								None then all subtypes are kept.
	@param PICT_FILE:		The path to a file which will be made into a pie
								chart of the subtypes.
	
	@returns				A list of ViralSeq objects
	"""
	
	output_list = []
	subtype_count = defaultdict(int)
	
	
	with open(INPUT_FILE) as handle:
		for this_seq in SeqIO.parse(handle, 'fasta'):
			this_virus = PyVirus.ViralSeq(this_seq, None)
			
			########DEBUG VALUE!!!!!!!!!!!!!!!!!!
			this_virus.tested_subtype = 'DEBUG'
			
			
			logging.debug('read sequence:' + this_virus.seq_name)
			subtype_count[this_virus.tested_subtype] += 1
			if WANTED_SUBS == None:
				output_list.append(this_virus)
			elif this_virus.tested_subtype in WANTED_SUBS:
				output_list.append(this_virus)
	
	logging.debug('Making Figure:' + PICT_FILE)
	pylab.pie(subtype_count.values(), labels = subtype_count.keys())
	pylab.savefig(PICT_FILE)
	
	
	return output_list
	
def ParseOptions(INPUT_STRING):
	"""
	Parses the input options
	
	"""
	
	
	
	parser = optparse.OptionParser()
	
	parser.add_option('-d', '--direc',
						dest = 'cache_dir',
						type = 'string',
						help = 'A working directory for file storage',
						default = None)
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
						help = 'Log File Name',
						default = None)
	parser.add_option('--config_file',
						dest = 'config_name',
						type = 'string',
						help = 'Config File Name, Overrides commandline input',
						default = None)
	
	#take input sequences
	seq_group = optparse.OptionGroup(parser, 'Seqeunce Inputs')
	
	seq_group.add_option('-i', '--input_fasta',
						dest = 'input_fasta',
						type = 'string',
						action = 'append',
						help = 'Whole Genome Sequences to visualize',
						default = None)
	seq_group.add_option('--ref_direc',
						dest = 'ref_direc',
						type = 'string',
						default = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\',
						help = 'Directory of Reference Sequences')
	seq_group.add_option('-s','--subtype',
						dest = 'wanted_subs',
						action = 'append',
			help = 'The desired subtypes in the output: can be used multiple times',
						default = None)
	
	parser.add_option_group(seq_group)
	
	pat_group = optparse.OptionGroup(parser, 'Patient Options')
	
	pat_group.add_option('--pat_direc',
					dest = 'pat_direc',
					default = None,
					type = 'string',
					help = 'The full path to a directory of patient sequences')
	pat_group.add_option('--sd_method',
					dest = 'sd_method',
					default = False,
					action = 'store_true',
					help = 'Use the SD method for determining responders')
					
	
	
	
	
	#options for types of outputs
	feat_group = optparse.OptionGroup(parser, 'Feature Outputs')
	feat_group.add_option('--ELM',
							dest = 'elm_flag',
							action = 'store_true',
							default = False,
							help = 'Toggle ELM Features')
	feat_group.add_option('--TF',
							dest = 'tf_flag',
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
							help = 'A file of known Homology Islands',
							default = None)
	
	parser.add_option_group(feat_group)
	
	#options for defining externa dependancies
	dep_group = optparse.OptionGroup(parser, 'Dependency Locations')
	
	dep_group.add_option('--blast_dir',
						dest = 'blast_dir',
						default = os.environ['BLASTPATH'],
						type = 'string',
						help = 'Location of Local Blast binaries')
	dep_group.add_option('--hybrid_dir',
						dest = 'hybrid_dir',
						default = os.environ['HYBRIDPATH'],
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
	
	fig_group = optparse.OptionGroup(parser, 'Figure Locations')
	
	fig_group.add_option('--out_config',
							dest = 'out_config',
							default = None,
							type = 'string',
							help = 'Location to output a lsit of figure locations')
	
	fig_group.add_option('--subtype_pie',
							dest = 'sub_figure',
							default = None,
							type = 'string',
							help = 'Location to output subtype Pie Chart')
	
	fig_group.add_option('--alignment_fig',
							dest = 'align_fig',
							default = None,
							type = 'string',
							help = 'Location to output Alignment Figure')
	
	fig_group.add_option('--clin_fig',
							dest = 'clin_fig',
							default = None,
							type = 'string',
					help = 'Location to output the Clinical Timecourse figure')
	
	parser.add_option_group(fig_group)
	
	out_feat_group = optparse.OptionGroup(parser, 'Feature Output')
	
	out_feat_group.add_option('--base_name',
								dest = 'base_name',
								default = 'feature',
								type = 'string',
								help = 'The base name to save feature-lists')
	
	parser.add_option_group(out_feat_group)
	
	
	(options, args) = parser.parse_args(INPUT_STRING)
	
	
	if options.config_name != None:
		this_config = ConfigParser.ConfigParser()
		this_config.read(options.config_name)
		config_dict = this_config.defaults()
		
		for this_key in config_dict:
			logging.warning(this_key + str(config_dict[this_key]))
			parser.set_default(this_key, config_dict[this_key])
		
		(options, args) = parser.parse_args(INPUT_STRING)
	
	logging.warning(str(options.num_threads))
	
	return options, args


	
class VisualController():
	def __init__(self, OPTIONS, ARGS):
		"""
		Initializes the controller class using the output created by optparse.
		"""
		
		self.options = OPTIONS
		self.args = ARGS
		self.outconfig = ConfigParser.ConfigParser()
		
		self.test_seqs = None
		self.pat_base = None
		
		
		self.LoadData()
		self.ProcessOptions()
	
	def LoadData(self):
		"""
		Loads the data into memory based on the options provided.
		"""
		
		logging.warning('Reading sequences')
		self.test_seqs = []
		if type(self.options.input_fasta) == types.ListType:
			for this_file in self.options.input_fasta:
				self.test_seqs += Fasta2ViralSeq(this_file, 
												self.options.wanted_subs,
												self.options.sub_figure)
		elif type(self.options.input_fasta) == types.StringType:
			self.test_seqs = Fasta2ViralSeq(self.options.input_fasta, 
												self.options.wanted_subs,
												self.options.sub_figure)
		
		logging.warning('%(num)d sequences read' % {'num': len(self.test_seqs)})
		#log the output
		logging.debug('Logging the output')
		self.outconfig.set('DEFAULT', 'sub_figure', self.options.sub_figure)
		self.outconfig.set('DEFAULT', 'wanted_seqs', len(self.test_seqs))
		
		if self.options.pat_direc != None:
			logging.warning('Reading Patient Sequences')
			self.pat_base = PatUtils.PatBase()
			self.pat_base.ReadDirec(self.options.pat_direc)
			logging.warning('Done reading Patient Sequences')
	
	
	def ProcessOptions(self):
		"""
		Processes the options provided
		"""
		
		
		#process reference sequences
		
		logging.warning('Creating MappingBase')
		self.mapping_base = HIVDatabase.MappingBase(self.options.ref_direc,
												self.options.scratch_dir, 
												'SEQ_FILE')
		
		logging.warning('Building Fasta Database')
		self.mapping_base.BuildRefBase()
		
		logging.warning('Adding sequences to the Shelf')
		if len(self.test_seqs) != 0:
			self.mapping_base.AddtoShelf(self.test_seqs)
		
		if self.pat_base != None:
			logging.warning('Adding Patient Sequences to the shelf')
			self.mapping_base.AddtoShelf(self.pat_base.values())
		if self.options.clin_fig != None:
			self.pat_base.PatTimeCourseFig(self.options.clin_fig)
			self.outconfig.set('DEFAULT', 'ClinicalTimcourseFigure', 
								self.options.clin_fig)
			
		
		
		
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
			# for this_calib in calib_dict.keys()[30:]:
				# junk = calib_dict.pop(this_calib)
		
		logging.warning('Beginning Annotation')
		
		if self.options.mirna_flag:
			ref_base.ref_seqs[0].HumanMiRNAsite(calib_dict)
		if self.options.elm_flag:
			ref_base.ref_seqs[0].FindELMs(elm_dict)
		if self.options.tf_flag:
			ref_base.ref_seqs[0].FindTFSites()
		if self.options.hom_flag:
			ref_base.ref_seqs[0].FindHomIslands(ref_base.ref_seqs[0])
		
		
		self.mapping_base.AnnotateBase(calib_dict, elm_dict, self.options.tf_flag, 
									self.options.hom_flag, True, 
									ref_base.ref_seqs[0].seq_name)
		
		
		resp_name = self.options.out_direc + self.options.base_name + '_resp.txt'
		nonresp_name = self.options.out_direc + self.options.base_name + '_nonresp.txt'
		annot_name = self.options.out_direc + self.options.base_name + '_description.txt'
		
		
		#color_list = map(lambda x: x.DetermineResponder('WND'), self.pat_base.values())
		
		t_resp_count = 0
		t_non_resp = 0
		resp = []
		non_resp = []
		for this_val in self.pat_base.values():
			if this_val.DetermineResponder('WND'):
				t_resp_count += 1
				resp.append(this_val)
			else:
				t_non_resp += 1
				non_resp.append(this_val)
		logging.debug('Test Found %(r)d responders and %(nr)d non-resp' % \
						{'r':t_resp_count, 'nr':t_non_resp})
		
		
		with open(resp_name, mode = 'w') as handle:
			for this_pat in resp:
				this_pat.WriteFeatures(ref_base.ref_seqs[0].feature_annot, handle)
		
		with open(nonresp_name, mode = 'w') as handle:
			for this_pat in non_resp:
				this_pat.WriteFeatures(ref_base.ref_seqs[0].feature_annot, handle)
		
		with open(annot_name, mode = 'w') as handle:
			for this_feat in ref_base.ref_seqs[0].feature_annot:
				handle.write(str(this_feat) + '\n')
		
		
		#filter_fun = lambda x: x.CheckRange('HumanMiRNA', None)
		filter_fun = None
		
		resp_fun = lambda x:x.DetermineResponder('SD')
		
		gene_fig = self.mapping_base.MakeMultiDiagram(ref_base.ref_seqs[0].seq_name,
														self.options.wanted_subs,
														ANCHOR_FILT = filter_fun,
														DISPLAY_GROUPING = resp_fun)
		
		if self.options.align_fig != None:
			gene_fig.draw(format = 'linear', fragments = 1)
			gene_fig.write(self.options.align_fig, 'PDF')
		
		self.outconfig.set('DEFAULT', 'AlignmentFigure', self.options.align_fig)
		
		

	def FinalCleanup(self):
		"""
		Cleans up the data afterwards and writes a configuration file 
		explaining the output.
		"""
		
		with open(self.options.out_config, mode = 'w') as handle:
			self.outconfig.write(handle)



if __name__ == '__main__':
	
	console_fmt = '%(asctime)-7s %(message)s'
	logging.basicConfig(level = logging.DEBUG, format = console_fmt)
	options, args = ParseOptions(sys.argv[1:])
	
	file_fmt = '%(levelname)s \t %(threadName)s \t %(funcName)s \t %(asctime)s \t %(message)s'
	
	
	
	#logging.addLevelName('CALC_RESULT', 15)
	
	#if options.verbose:
	
	#else:
	#	logging.basicConfig(level = logging.WARNING, 
	#						format = console_fmt)
	
	# console = logging.StreamHandler()
	# console.setLevel(level = logging.DEBUG)
	# cons_formatter = logging.Formatter(console_fmt)
	# console.setFormatter(cons_formatter)
	# logging.getLogger('').addHandler(console)
	
	
	
	logging.debug('correct console adding')
	if options.log_name != None:
		file_name = options.out_direc + options.log_name
		file_log = logging.FileHandler(file_name, mode = 'w')
		formatter = logging.Formatter(file_fmt)
		file_log.setLevel(level = logging.DEBUG)
		file_log.setFormatter(formatter)
		logging.getLogger('').addHandler(file_log)
	
	logging.info('Finished Input Parsing')
	
	
	this_visual = VisualController(options, args)
	
	logging.info('Finished Pre-Processing')
	
	
	this_visual.run()
	
	logging.info('Finished Running')
	
	this_visual.FinalCleanup()
	
	
	#start doing stuff with input options
	
	
	
	
	
	
	
	
	
	


