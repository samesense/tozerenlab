from __future__ import with_statement
import itertools as IT
import HIVGenoTypeChecker as GenoTyping
import subprocess
import tempfile
import threading
import time
import os
import re
import types
import Queue
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein, generic_nucleotide
from Bio.Blast import NCBIXML, NCBIStandalone
import ctypes


def fLock(FILENAME):
	print 'in fLock'
	h = ctypes.cdll.msvcrt._sopen(FILENAME, 0, 0x10, 0)
	while h == -1:
		time.sleep(0.5)
		h = ctypes.cdll.msvcrt._sopen(FILENAME, 0, 0x10, 0)
	ctypes.cdll.msvcrt._close(h)
	print 'leaving fLock'


class Gene():
	def __init__(self, NAME, PRODUCT, NT_SEQ, AA_SEQ, START, END):
		self.nt_seq = NT_SEQ
		self.aa_seq = AA_SEQ
		self.start = START
		self.end = END
		self.name = NAME
		self.product = PRODUCT

	def __str__(self):
		temp_str = '(Gene Class '
		temp_str += 'Name: ' + self.name
		temp_str += ', Start: ' + str(self.start) + ')'
		return temp_str
	def __hash__(self):
		return hash(self.aa_seq)


class ViralSeq():
	def __init__(self, TEST_SEQ, SEQ_NAME):
		#if SEQ_NAME is None then assume that TEST_SEQ is a SeqRecord
		if SEQ_NAME == None:
			self.my_sequence = TEST_SEQ.seq.tostring()
			self.seq_name = TEST_SEQ.id
		else:
			self.my_sequence = TEST_SEQ
			self.seq_name = SEQ_NAME
		self.tested_subtype = None
		self.annotation = None

	def __hash__(self):
		return hash(self.my_sequence + self.seq_name)

	def GenomeToBioPython(self):
		"""
		GenomeToBioPython
				Returns a generic SeqRecord in the biopython format for
		the whole genome sequence.
		"""
		return SeqRecord(Seq(self.my_sequence, generic_nucleotide),
						 id=self.seq_name)

	def AnnotToBioPython(self):
		"""
		GenomeToBioPython
				Returns a list of generic SeqRecords in the biopython
		format for each protein in self.annotation.
		"""
		final_list = []
		for this_gene in self.annotation:
			this_record = SeqRecord(Seq(self.annotation[this_gene].aa_seq,
										generic_protein),
										id = self.seq_name + ':' + this_gene)
			final_list.append(this_record)
		return final_list


	def DetSubtype(self):
		"""
		DetSubtype
				Determines the subtype of the sequence using the
		HIVSubtype module.
		"""

		self.tested_subtype = GenoTyping.GetSimple(self.my_sequence)

	def TranslateAll(self, REFBASE):
		"""
		TranslateAll
				Uses a BLASTx query to determine the which ORFs
		correspond to HIV-1 proteins.  It then stores the data in a
		dictionary keyed by "genename".
		"""
		blast_record = REFBASE.BLASTx(self.GenomeToBioPython())
		gene_dict = {}
		reg = re.compile('.*?:(\w{3}).*')
		for this_check in blast_record.alignments:
			this_gene = str(reg.match(this_check.title).groups()[0])
			if not(gene_dict.has_key(this_gene)):
				aa_seq = str(this_check.hsps[0].query)
				start_pos = this_check.hsps[0].query_start
				end_pos = start_pos + this_check.length*3
				nt_seq = self.my_sequence[start_pos:end_pos]
				gene = Gene(this_gene, this_gene.upper(),
							nt_seq, aa_seq, start_pos, end_pos)
				gene_dict[this_gene] = gene
		self.annotation = gene_dict

class BkgSeq(ViralSeq):
	"""
	Used for Background Sequences
	"""


class PatSeq(ViralSeq):
	def __init__(self, TEST_SEQ, SEQ_NAME):
		self.my_sequence = TEST_SEQ
		self.pat_data = None
		self.tested_subtype = None
		self.seq_name = SEQ_NAME


class RefSeq(ViralSeq):
	def __init__(self, FILENAME):
		self.my_sequence = None
		self.known_subtype = None
		self.seq_name = None
		self.annotation = None

		self.ParseGenbank(FILENAME)


	def ParseGenbank(self, FILENAME):
		with open(FILENAME, mode='r') as handle:
			this_record = SeqIO.parse(handle, 'genbank').next()

		self.my_sequence = this_record.seq.tostring()
		self.seq_name = this_record.id
		

		self.annotation = {}
		for feat in this_record.features:
			if feat.qualifiers.has_key('gene') & feat.qualifiers.has_key('translation'):
				gene_name = feat.qualifiers['gene'][0]
				if '-' in gene_name:
					continue
				start_pos = feat.location._start.position
				end_pos = feat.location._end.position
				trans_data = feat.qualifiers['translation'][0]
				nt_seq = self.my_sequence[start_pos:end_pos]
				
				prod_name = feat.qualifiers['product'][0]
				
				this_gene = Gene(gene_name, prod_name, nt_seq,
								 trans_data, start_pos, end_pos)

				self.annotation[feat.qualifiers['gene'][0]] = this_gene



class RefBase():
	def __init__(self, SOURCE_DIR, DEST_DIR, BUILD = False,
				 BLAST_DIR = 'C:\\local_blast\\'):


		self.ref_seqs = []
		for this_file in filter(lambda x: x[-3:] == '.gb',
								os.listdir(SOURCE_DIR)):
				
			self.ref_seqs.append(RefSeq(SOURCE_DIR + this_file))
		
		self.nt_name = DEST_DIR + 'ref_nt.fasta'
		self.aa_name = DEST_DIR + 'ref_aa.fasta'
		self.dir_name = DEST_DIR
		self.blast_dir = BLAST_DIR

		if BUILD:
			self.BuildDatabase()

	def BuildDatabase(self):
		"""
		BuildDatabase
				Uses the RefSeqs provided to create a protein and nucleotide
		blast database.
		"""

		this_iter = IT.imap(RefSeq.GenomeToBioPython, iter(self.ref_seqs))
		with open(self.nt_name, mode='w') as handle:
			SeqIO.write(this_iter, handle, "fasta")

		with open(self.aa_name, mode='w') as handle:
			for this_data in self.ref_seqs:
				this_iter = this_data.AnnotToBioPython()
				SeqIO.write(this_iter, handle, "fasta")

		command = self.MakeBLASTCommand('formatdb',
										self.nt_name,
										EXTRA_FLAGS = ' -p F -o T')
		
		ProcessVar = subprocess.Popen(command, shell = True)
		ProcessVar.wait()

		command = self.MakeBLASTCommand('formatdb',
										self.aa_name,
										EXTRA_FLAGS = ' -p T -o T')

		ProcessVar = subprocess.Popen(command, shell = True)
		ProcessVar.wait()

	def BLASTx(self, INPUT_SEQ_RECORD, NUM_THREADS = None):
		"""
		BLASTx
				Performs a BLASTx query on the database to determine
		the possible translations of the INPUT_GENOME
		"""
		if NUM_THREADS == None:
			b_controler = BLASTController(INPUT_SEQ_RECORD, 'BLASTx', NUM_SEQS = 1, 
										BLAST_DIR = self.blast_dir,
										SCRATCH_DIR = self.dir_name)
			return b_controler.GetSingleResult()
		else:
			b_controler = BLASTController(INPUT_SEQ_RECORD, 'BLASTx', 
										BLAST_DIR = self.blast_dir,
										SCRATCH_DIR = self.dir_name, 
										NUM_THREADS = NUM_THREADS)
			b_controler.start()
			return b_controler.ResultGen()

	def BLASTn(self, INPUT_SEQ_RECORD, NUM_THREADS = None):
		"""
		BLASTn
				Performs a BLASTn query to align the provided genome to
		the reference genomes.
		"""
		if NUM_THREADS == None:
			b_controler = BLASTController(INPUT_SEQ_RECORD, 'BLASTn', NUM_SEQS = 1, 
										BLAST_DIR = self.blast_dir,
										SCRATCH_DIR = self.dir_name)
			return b_controler.GetSingleResult()
		else:
			b_controler = BLASTController(INPUT_SEQ_RECORD, 'BLASTn', 
										BLAST_DIR = self.blast_dir,
										SCRATCH_DIR = self.dir_name, 
										NUM_THREADS = NUM_THREADS)
			b_controler.start()
			return b_controler.ResultGen()
			
	def MakeBLASTCommand(self, BLAST_TYPE, INPUT_FILE, EXTRA_FLAGS = None):
		"""
		MakeBLASTCommand
				Creates a command that can be passed to Popen for
		executing blast commands ... supports: formatdb, blastn,
		blastx.
		"""
		if BLAST_TYPE == 'formatdb':
			command = self.blast_dir + 'bin\\formatdb.exe'

		command += ' -i ' + INPUT_FILE
		
		if EXTRA_FLAGS != None:
			command +=  EXTRA_FLAGS

		return command



class BLASTContext:
        """
        Provides a context where the temp_files are removed upon exit.
        """
	def __init__(self, BASE_DIR, NUM_FILES, PREFIX = 'temp_'):
		self.base_dir = BASE_DIR
		self.file_names = []
		for i in xrange(NUM_FILES):
			 newNamedFile = tempfile.NamedTemporaryFile(prefix = PREFIX,
											   dir = BASE_DIR)
			 self.file_names.append(newNamedFile.name)
			 newNamedFile.close()
			
	def __enter__(self):
		return self.file_names

	def __exit__(self, TYPE_IN, VALUE_IN, TRACEBACK_IN):
		for this_file in self.file_names:
			try:
				os.remove(this_file)
			except:
				continue

class BLASTController:

	def __init__(self, SEQ_ITER, BLAST_TYPE, NUM_SEQS = None,
				 BLAST_DIR = "C:\\local_blast\\",
				 SCRATCH_DIR = "C:\\pyscratch\\",
				 NUM_THREADS = 4, AA_BASE = 'ref_aa.fasta',
				 NT_BASE = 'ref_nt.fasta'):
		if BLAST_TYPE == 'BLASTn':
			self.blast_type = 'blastn'
		elif BLAST_TYPE == 'BLASTx':
			self.blast_type = 'blastx'
		else:
			raise KeyError, 'Unrecognized BLAST type'
		if NUM_SEQS == 1:
			self.seqs = iter([SEQ_ITER])
		else:
			self.seqs = SEQ_ITER
			
		self.thread_seph = threading.Semaphore(NUM_THREADS)
		self.file_seph = threading.Semaphore(20)
		self.seq_file_names = Queue.Queue()
		self.res_file_names = Queue.Queue()
		self.scratch_dir = SCRATCH_DIR
		self.blast_dir = BLAST_DIR
		self.aa_name = SCRATCH_DIR + AA_BASE
		self.nt_name = SCRATCH_DIR + NT_BASE

	def start(self):
		"""
		start
			Starts the multithreaded BLASTing
		"""
		cont_thread = threading.Thread(name = 'BLAST-CONT',
						target = self.runBLASTs)
		cont_thread.start()
					
	def runBLASTs(self):
		"""
		runBLASTs
			Controls the starting and stopping of the workers
		"""
		this_seq = self.seqs.next()
		these_files = self.MakeFiles()
		with open(these_files[0], mode = 'w') as seq_handle:
			SeqIO.write([this_seq],
						seq_handle, 'fasta')
		#wait until there is an availble thread
		self.thread_seph.acquire()
		this_thread = self.BLASTWorker(self.blast_type, these_files[0],
						these_files[1], self.blast_dir, self.nt_name,
						self.aa_name, self.thread_seph)
		this_thread.start()
		this_thread.join()
		self.seq_file_names.put(these_files[0])
		self.res_file_names.put(these_files[1])
	
		
		
		for this_seq in self.seqs:
			#autimatically limits the number of files possible
			these_files = self.MakeFiles()
			with open(these_files[0], mode = 'w') as seq_handle:
				SeqIO.write([this_seq],
							seq_handle, 'fasta')
			#wait until there is an availble thread
			self.thread_seph.acquire()
			this_thread = self.BLASTWorker(self.blast_type, these_files[0],
							these_files[1], self.blast_dir, self.nt_name,
							self.aa_name, self.thread_seph)
			this_thread.start()
			self.seq_file_names.put(these_files[0])
			self.res_file_names.put(these_files[1])
				
		self.seq_file_names.put(StopIteration)
		self.res_file_names.put(StopIteration)
		
	def ResultGen(self):
		"""
		ResultGen
			A generator function which returns the blast results in the 
			order in which they were provided.
		"""
		while True:
			this_seq_file = self.seq_file_names.get()
			this_res_file = self.res_file_names.get()
			if this_seq_file == StopIteration:
				break
			
			fLock(this_res_file)
			
			with open(this_res_file, mode = 'r') as handle:
				blast_rec = NCBIXML.parse(handle).next()
				
			fLock(this_res_file)
			os.remove(this_res_file)
			self.file_seph.release()
			
			fLock(this_seq_file)
			os.remove(this_seq_file)
			self.file_seph.release()
			
			
			yield blast_rec
				
			
	def MakeFiles(self):
		"""
		MakeFiles
				Creates 2 temporary files (one for sequence and one
				for blast result).  Uses the self.file_seph to ensure
				that only a limited number of files are created at a
				time.
		"""
		file_names = []
		for i in xrange(2):
			self.file_seph.acquire()
			newNamedFile = tempfile.NamedTemporaryFile(prefix = 'temp_', 
														dir = self.scratch_dir)
			file_names.append(newNamedFile.name)
			newNamedFile.close()
		return file_names


	def GetSingleResult(self):
		"""
		GetSingleResult
			A shortcut for getting a single BLAST result without any threading
			nonsense.
		"""
		
		these_files = self.MakeFiles()
		with open(these_files[0], mode = 'w') as seq_handle:
				SeqIO.write([self.seqs.next()],
							seq_handle, 'fasta')
		
		this_worker = self.BLASTWorker(self.blast_type, these_files[0],
							these_files[1], self.blast_dir, self.nt_name,
							self.aa_name, self.thread_seph)
		
		command = this_worker.MakeBLASTCommand()
		ProcessVar = subprocess.Popen(command, shell = True)
		ProcessVar.wait()
		
		with open(these_files[1], mode = 'r') as handle:
			blast_rec = NCBIXML.parse(handle).next()
		
		
		os.remove(these_files[0])
		os.remove(these_files[1])
		
		return blast_rec
					
		
	class BLASTWorker(threading.Thread):
		def __init__(self, BLAST_TYPE, SEQ_FILE, OUTPUT_FILE,
					 BLAST_DIR, NT_NAME, AA_NAME, THREAD_SEPH):
			if BLAST_TYPE == 'blastn':
				self.blast_type = 'blastn'
			elif BLAST_TYPE == 'blastx':
				self.blast_type = 'blastx'
			else:
				raise KeyError, 'Unrecognized BLAST type'
				
			threading.Thread.__init__(self)
			self.seq_file = SEQ_FILE
			self.output_file = OUTPUT_FILE
			self.blast_dir = BLAST_DIR
			self.aa_name = AA_NAME
			self.nt_name = NT_NAME
			self.thread_seph = THREAD_SEPH


		def run(self):
			"""
			run
				Performs the BLAST query and waits for
				completion.
			"""
			
			command = self.MakeBLASTCommand()
			ProcessVar = subprocess.Popen(command, shell = True)
			ProcessVar.wait()
			self.thread_seph.release()
			
				   
						   
		def MakeBLASTCommand(self):
			"""
			MakeBLASTCommand
					Creates a command that can be passed to Popen for
			executing blast commands ... supports: formatdb, blastn,
			blastx.
			"""

			if self.blast_type == 'blastx':
				command = self.blast_dir + 'bin\\blastall.exe'
				command += ' -p blastx -m 7'
				command += ' -d ' + self.aa_name
				command += ' -o ' + self.output_file
					
			if self.blast_type == 'blastn':
				command = self.blast_dir + 'bin\\blastall.exe'
				command += ' -p blastn -m 7'
				command += ' -d ' + self.nt_name
				command += ' -g F'
				command += ' -o ' + self.output_file

			command += ' -i ' + self.seq_file
			

			return command




        
        


    
