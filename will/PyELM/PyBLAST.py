from __future__ import with_statement
import tempfile
import threading
import subprocess
from Bio import SeqIO
from Bio.Blast import NCBIXML, NCBIStandalone
import Queue
import ctypes
import os
import time
import logging

logger = logging.getLogger('')

def fLock(FILENAME):
	logger.debug('Getting fLock:' + FILENAME)
	h = ctypes.cdll.msvcrt._sopen(FILENAME, 0, 0x10, 0)
	while h == -1:
		time.sleep(0.5)
		h = ctypes.cdll.msvcrt._sopen(FILENAME, 0, 0x10, 0)
	ctypes.cdll.msvcrt._close(h)
	logger.debug('Releasing fLock:' + FILENAME)



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
		oounter = 0
		this_seq = self.seqs.next()
		oounter += 1
		these_files = self.MakeFiles(oounter)
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
			oounter += 1
			#autimatically limits the number of files possible
			these_files = self.MakeFiles(oounter)
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
				
			
	def MakeFiles(self, NUM):
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
			pref = 'temp_' + str(NUM) + '_'
			newNamedFile = tempfile.NamedTemporaryFile(prefix = pref, 
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
		
		these_files = self.MakeFiles(1)
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
		#time.sleep(0.5)
		fLock(these_files[0])
		os.remove(these_files[0])
		fLock(these_files[1])
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
			logger.debug('BLAST command sent')
			ProcessVar.wait()
			self.thread_seph.release()
			logger.debug('BLAST finished')
			
				   
						   
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
					
			if self.blast_type == 'blastn':
				command = self.blast_dir + 'bin\\blastall.exe'
				command += ' -p blastn -m 7'
				command += ' -d ' + self.nt_name
			
			command += ' -o ' + self.output_file
			command += ' -i ' + self.seq_file
			

			return command
