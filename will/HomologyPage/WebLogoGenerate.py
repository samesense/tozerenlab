from __future__ import with_statement
import weblogolib
import tempfile
import os


class FileContext:
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

def GenerateSeqLogo(SEQS, OUT_FILE):
	"""
	Generates a SeqLogo based on the input sequences and saves the image into 
	OUT_FILE
	
	@param SEQ:			A list of sequences.
	@param OUT_FILE:	The path to the output image location.
	"""
	with FileContext(os.curdir, 1, PREFIX = 'img') as f_name:
		with open(f_name[0], mode = 'w') as handle:
			for this_seq in SEQS:
				handle.write('>test_seq\n' + this_seq + '\n')
		with open(f_name[0]) as handle:
			read_seqs = weblogolib.read_seq_data(handle)
	data = weblogolib.LogoData.from_seqs(read_seqs)
	options = weblogolib.LogoOptions(unit_name = 'probability')
	
	format = weblogolib.LogoFormat(data, options)
	with open(OUT_FILE, mode = 'w') as handle:
		weblogolib.eps_formatter(data, format, handle)
	