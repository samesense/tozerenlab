from __future__ import with_statement
import re
import PyMozilla
import os
import shelve
import string
from Bio import SeqIO
import threading
import logging
import itertools
import copy
import optparse



class Annot():
	def __init__(self, NAME, START_POS, END_POS, TYPE):
		self.name = NAME
		self.start = START_POS
		self.end = END_POS
		self.type = TYPE
	
	def __eq__(self, VAL):
		return self.name == VAL.name
		
	def __str__(self):
		this_str = self.name + '\t'
		this_str += str(self.start) + '\t'
		this_str += str(self.end)
		
		return this_str


class Protein():
	def __init__(self, REF_PROT, AA_SEQ):
		self.aa_seq = AA_SEQ
		self.ref_prot = REF_PROT
		self.annot = []
		#self.neighbors = []
		
	def HasAnnot(self, ANNOT_NAME):
		"""
		Checks to see if this Protein has requested ELM
		"""
		return Annot(ANNOT_NAME, None, None, None) in self.annot
		
		
	def AnnotELMs(self, CACHE_DIR, FORCE_DL = 0):
		"""
		Makes a call to CheckELM and ReadData to annotate the ELMs on this
		protein.  Saves the downloaded .html file to CACHE_DIR to speed up 
		future retrievals.
		FORCE_DL == 0
			Will only process the cached data and NEVER make a server request.
		FORCE_DL == 1
			Will process cached data if present and will make a server request
			if the data is absent in the cache.
		FORCE_DL == 2
			Will ALWAYS make a server request and will overwrite the cached 
			data.
		
		"""
		
		if FORCE_DL == 0:
			self.annot = ReadData(CACHE_DIR + self.ref_prot + '.html')
			logging.info('File Read\t ' + self.ref_prot)
			#this_string = string.join(self.elm_annot())
			#logging.info()
			return
		
		
		file_present = (self.ref_prot + '.html') in os.listdir(CACHE_DIR)
		want_retrieve = FORCE_DL >= 1
		force_retrive = FORCE_DL >= 2
		
		if force_retrive | (want_retrieve and not file_present):
			try:
				CheckELM(self.aa_seq, CACHE_DIR + self.ref_prot + '.html')
			except:
				logging.warning('Could not retrive\t ' + self.ref_prot)
		
		self.annot = ReadData(CACHE_DIR + self.ref_prot + '.html')
		logging.info('File Read\t ' + self.ref_prot)

		
def ReadData(INPUT_FILENAME):
	other_re = re.compile('<TR><TD valign=top><A name=".*?:in">.*?</TD></TR>', re.S)
	name_re = re.compile('<A name="(.*?):in"')
	spot_re = re.compile('<td>\d{1,}-\d{1,}</td>')
	num_re = re.compile('\d{1,}')
	
	
	elm_vec = []
	try:
		with open(INPUT_FILENAME) as handle:
			html_data = handle.read()
	except:
		name = string.split(INPUT_FILENAME, os.sep)[-1]
		logging.warning('Missing file\t' + name)
		return elm_vec
	
	
	for this_check in other_re.findall(html_data):
		elm_name = name_re.search(this_check).groups()[0]
		for this_num in spot_re.findall(this_check):
			nums = num_re.findall(this_num)
			elm_vec.append(Annot(elm_name, int(nums[0]), int(nums[1]), 'ELM'))
	
	return elm_vec
	
	
	
def Fasta2Shelve(FASTA_F_NAME, IND_FUN = None):
	"""
	Converts a fasta file into a shelve that is indexed by the fasta-header.
	If a IND_FUN is provided then the output of that function will be used as
	this index key.
	"""
	
	fasta_shelve = {}
	
	with open(FASTA_F_NAME) as fasta_handle:
		seq_iter = SeqIO.parse(fasta_handle, 'fasta')
		
		for this_seq in seq_iter:
			this_id = this_seq.id
			if IND_FUN != None:
				this_id = IND_FUN(this_id)
			this_letters = this_seq.seq.tostring()
			fasta_shelve[this_id] = this_letters
			logging.info('Shelved\t ' + this_id)
			
	return fasta_shelve

def CheckELM(INPUT_SEQ, OUTPUT_FILENAME):
	moz_emu = PyMozilla.MozillaEmulator(cacher=None)
	base_url = 'http://elm.eu.org/basicELM/'
	redirect_re = re.compile('<META HTTP-EQUIV="REFRESH" CONTENT="10; URL=(.*?);r=1">')
	input_data = [('sequence', INPUT_SEQ), ('userSpecies', '9606'), ('fun', 'Submit')]
	
	name = string.split(OUTPUT_FILENAME, os.sep)[-1]
	
	logging.info('Made Server Request\t ' + name)
	
	output = moz_emu.post_multipart(base_url+'cgimodel.py',input_data, 
									[], forbid_redirect=False)
	#time.sleep(10)
	new_url = base_url + redirect_re.search(output).groups()[0]
	sec_out = moz_emu.download(new_url)
	with open(OUTPUT_FILENAME, mode='w') as handle:
		handle.write(sec_out)
	
	logging.info('Data written\t ' + name)

def AnnotWorker(SEQ_ID, AA_SEQ, CACHE_DIR):
	"""
	The threaded worker function for annotation
	"""
	
	this_pro = Protein(SEQ_ID, AA_SEQ)
	this_pro.AnnotELMs(CACHE_DIR, FORCE_DL = 1)
	
	with elm_dict_lock:
		elm_dict[SEQ_ID] = this_pro
	
	logging.info('Logged\t ' + SEQ_ID)
	worker_seph.release()
	
	
	
if __name__ == '__main__':
	
	parser = optparse.OptionParser()
	parser.add_option('-i', '--input_file',
						dest = 'INPUT_FASTA_FILE',
						help = 'Input File Path')
	parser.add_option('-o', '--output_file',
						dest = 'OUTPUT_FILE',
						default = 'elm_results.txt',
						help = 'Output File Path')
	parser.add_option('-d', '--cache_dir',
						dest = 'CACHE_DIR',
						default = os.curdir,
						help = 'Directory to save cached results')
	parser.add_option('-v', '--verbose',
						dest = 'VERBOSITY',
						action = 'store_true',
						default = False,
						help = 'Be Verbose')
	parser.add_option('-f', '--force_dl',
						dest = 'FORCE_DOWNLOAD',
						default = False)
	parser.add_option('-l', '--log_file',
						dest = 'log_file',
						default = 'log_file.txt',
						help = 'Path to a log file')
	parser.add_option('-t', '--threads',
						dest = 'NUM_THREADS',
						default = 10,
						type = 'int',
						help = 'Number of Threads to use for the ELM server')
	
	options, remainder = parser.parse_args()
	log_file = options.log_file
	VERBOSITY = options.VERBOSITY
	CACHE_DIR = options.CACHE_DIR
	INPUT_FASTA_FILE = options.INPUT_FASTA_FILE
	OUTPUT_FILE = options.OUTPUT_FILE
	FORCE_DOWNLOAD = options.FORCE_DOWNLOAD
	NUM_THREADS = options.NUM_THREADS
	
	
	
	
	
	fmt='%(threadName)s \t %(funcName)s \t %(asctime)s \t %(message)s'
	logging.basicConfig(level=logging.DEBUG,
						filename=log_file, 
						filemode='w',
						 format = fmt)
	#logging.basicConfig(level=logging.INFO, format = fmt)
	console = logging.StreamHandler()
	if VERBOSITY:
		console.setLevel(logging.INFO)
	else:
		console.setLevel(logging.WARNING)
	
	formatter = logging.Formatter('%(name)-12s: %(message)s')
	console.setFormatter(formatter)
	logging.getLogger('').addHandler(console)
	
	seq_dict = Fasta2Shelve(INPUT_FASTA_FILE)
	
	all_threads = []
	
	worker_seph = threading.Semaphore(NUM_THREADS)
	
	elm_dict_lock = threading.Lock()
	elm_dict = {}
	
	
	
	for this_seq in seq_dict:
		worker_seph.acquire()
		this_thread = threading.Thread(target = AnnotWorker,
										args = (this_seq, seq_dict[this_seq], 
												CACHE_DIR),
										name = this_seq)
		all_threads.append(this_thread)
		this_thread.start()
		
	for this_thread in all_threads:
		this_thread.join()
		
	with open(OUTPUT_FILE, mode = 'w') as handle:
		for this_seq in elm_dict:
			for this_elm in elm_dict[this_seq].annot:
				handle.write(this_seq + '\t' + str(this_elm) + '\n')
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	