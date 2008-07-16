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
										self.nt_name, self.nt_name,
										EXTRA_FLAGS = ' -p F -o T')
		
		ProcessVar = subprocess.Popen(command, shell = True)
		ProcessVar.wait()

		command = self.MakeBLASTCommand('formatdb',
										self.aa_name, self.aa_name,
										EXTRA_FLAGS = ' -p T -o T')

		ProcessVar = subprocess.Popen(command, shell = True)
		ProcessVar.wait()

	def BLASTx(self, INPUT_GENOME):
		"""
		BLASTx
				Performs a BLASTx query on the database to determine
		the possible translations of the INPUT_GENOME
		"""
		with BLASTContext(self.dir_name, 2) as file_names:

			with open(file_names[0], mode = 'w') as handle:
				SeqIO.write([INPUT_GENOME],
							handle, 'fasta')

			command = self.MakeBLASTCommand('blastx', file_names[0],
											file_names[1])
			ProcessVar = subprocess.Popen(command)
			ProcessVar.wait()
		   
			with open(file_names[1], mode = 'r') as handle:
				this_blast = NCBIXML.parse(handle).next()

		return this_blast

	def BLASTn(self, INPUT_GENOME, DUMP_LOCATION = None):
		"""
		BLASTn
				Performs a BLASTn query to align the provided genome to
		the reference genomes.
		"""
		if type(INPUT_GENOME) != types.ListType:
			INPUT_GENOME = [INPUT_GENOME]

		with BLASTContext(self.dir_name, 2) as file_names:

			with open(file_names[0], mode = 'w') as handle:
				SeqIO.write(INPUT_GENOME,
							handle, 'fasta')

			command = self.MakeBLASTCommand('blastn', file_names[0],
											file_names[1],
											EXTRA_FLAGS = ' -B ' + str(len(INPUT_GENOME)))

			ProcessVar = subprocess.Popen(command, shell = True)
			ProcessVar.wait()

			with open(file_names[1], mode = 'r') as handle:
				blast_rec = NCBIXML.parse(handle).next()
		if DUMP_LOCATION == None:
			return blast_rec
		else:
			DUMP_LOCATION.put()
					

	def MultiThreadedBLAST(INPUT_SEQ_ITER, BLAST_TYPE = 'BLASTn',
						   NUM_THREADS = 4):
		"""
		MultiThreadedBLAST
				A controller for multithreading the BLAST calls.  This
		will leave each BLAST result in a temp_file and then remove them
		upon yield'ing the result.
		"""
		
		if BLAST_TYPE == 'BLASTn':
			worker_fun = lambda x: self.BLASTn(x)
		elif BLAST_TYPE == 'BLASTx':
			worker_fun = lambda x: self.BLASTx(x)
		else:
			raise KeyError, 'Unrecognized BLAST type'
		this_seph = threading.Semaphore(NUM_THREADS)
		this_ans = Queue.Queue(0)

		for this_seq in INPUT_SEQ_ITER:
			with this_seph:
				this_thread = threading.Thread(target = worker_fun,
											   args = (this_seq))
				this_thread.start()
	

	def MakeBLASTCommand(self, BLAST_TYPE, INPUT_FILE, OUTPUT_FILE,
						 EXTRA_FLAGS = None):
		"""
		MakeBLASTCommand
				Creates a command that can be passed to Popen for
		executing blast commands ... supports: formatdb, blastn,
		blastx.
		"""
		if BLAST_TYPE == 'formatdb':
			command = self.blast_dir + 'bin\\formatdb.exe'

		if BLAST_TYPE == 'blastx':
			command = self.blast_dir + 'bin\\blastall.exe'
			command += ' -p blastx -m 7'
			command += ' -d ' + self.aa_name
			command += ' -o ' + OUTPUT_FILE
			
		if BLAST_TYPE == 'blastn':
			command = self.blast_dir + 'bin\\blastall.exe'
			command += ' -p blastn -m 7'
			command += ' -d ' + self.nt_name
			command += ' -g F'
			command += ' -o ' + OUTPUT_FILE

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
				 BLAST_DIR = "C:\\local_blast\\bin\\",
				 SCRATCH_DIR = "C:\\pythonscrath\\",
				 NUM_THREADS = 4):
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
			
		self.thread_seph = threading.Sephamore(NUM_THREADS)
		self.file_seph = threading.Sephamore(min(50, NUM_SEQS))
		self.seq_file_names = Queue.Queue()
		self.res_file_names = Queue.Queue()
		self.scrath_dir = SCRATCH_DIR

	def start(self):
		"""
		start
			Starts the multithreaded BLASTing
		"""
		
		for this_seq in self.seqs:
			#autimatically limits the number of files possible
			these_files = self.MakeFiles()
			self.seq_file_names.put(these_files[0])
			self.res_file_names.put(these_files[1])
			#wait until there is an availble thread
			with self.thread_seph:
				self.BLASTWorker(self.blast_type, these_files[0],
								these_files[1], self.blast_dir)
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
			file_size = os.path.getsize(this_res_file)
			while file_size < 10:
				time.sleep(5)
				file_size = os.path.getsize(this_res_file)
			with open(file_names[1], mode = 'r') as handle:
				blast_rec = NCBIXML.parse(handle).next()
			
			self.file_seph.release()				
			os.remove(this_seq_file)
			self.file_seph.release()
			os.remove(this_res_file)
			
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
			self.file_seph.aquire()
			newNamedFile = tempfile.NamedTemporaryFile(prefix = PREFIX,
									   dir = self.scrath_dir)
			file_names.append(newNamedFile.name)
			newNamedFile.close()


	class BLASTWorker(threading.Thread):
		def __init__(self, BLAST_TYPE, SEQ_FILE, OUTPUT_FILE,
					 BLAST_DIR):
			if BLAST_TYPE == 'BLASTn':
				self.blast_type = 'blastn'
			elif BLAST_TYPE == 'BLASTx':
				self.blast_type = 'blastx'
			else:
				raise KeyError, 'Unrecognized BLAST type'
			self.seq_file = SEQ_FILE
			self.output_file = OUTPUT_FILE
			self.blast_dir = BLAST_DIR


		def run(self):
			"""
			run
				Performs the BLAST query and waits for
				completion.
			"""
			if BLAST_TYPE == 'BLASTn':
				blast_type = 'blastn'
			elif BLAST_TYPE == 'BLASTx':
				blast_type = 'blastx'
			else:
				raise KeyError, 'Unrecognized BLAST type'

			command = self.MakeBLASTCommand(selfl.blast_type,
											self.seq_file,
											self.output_file)
			ProcessVar = subprocess.Popen(command, shell = True)
			ProcessVar.wait()
			
				   
						   
		def MakeBLASTCommand(self, BLAST_TYPE, INPUT_FILE, OUTPUT_FILE,
							 EXTRA_FLAGS = None):
				"""
				MakeBLASTCommand
						Creates a command that can be passed to Popen for
				executing blast commands ... supports: formatdb, blastn,
				blastx.
				"""

			if BLAST_TYPE == 'blastx':
				command = self.blast_dir + 'bin\\blastall.exe'
				command += ' -p blastx -m 7'
				command += ' -d ' + self.aa_name
				command += ' -o ' + OUTPUT_FILE
					
			if BLAST_TYPE == 'blastn':
				command = self.blast_dir + 'bin\\blastall.exe'
				command += ' -p blastn -m 7'
				command += ' -d ' + self.nt_name
				command += ' -g F'
				command += ' -o ' + OUTPUT_FILE

			command += ' -i ' + INPUT_FILE
			
			if EXTRA_FLAGS != None:
				command +=  EXTRA_FLAGS

			return command




        
        


    
