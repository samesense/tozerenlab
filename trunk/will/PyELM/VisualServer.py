from __future__ import with_statement
import string,cgi,time,os,logging,re, PyMozilla, tempfile, ConfigParser
import subprocess

import itertools as IT
from Bio import SeqIO
from BaseHTTPServer import BaseHTTPRequestHandler, HTTPServer

logging.basicConfig(level=logging.DEBUG,
					format='%(asctime)s %(levelname)s %(message)s')

MakeFname = lambda x: string.join(x, os.sep)
					
					
class MyHandler(BaseHTTPRequestHandler):

	def do_GET(self):
		logging.debug('Recieved GET request')
		try:
			if self.path.endswith('VisualPage.html'):
				logging.debug('Correct PageName')
				f = open(DIREC + 'VisualPage.html') #self.path has /test.html
#note that this potentially makes every file on your computer readable by the internet

				self.send_response(200)
				self.send_header('Content-type', 'text/html')
				self.end_headers()
				self.wfile.write(f.read())
				f.close()
				return
			elif self.path.endswith('_resultsPage.html'):
				dir_name = re.findall('.*?(tmp.*?)_resultsPage.html', 
											self.path)
				
				f = open(os.environ['PYTHONSCRATCH'] + dir_name + 'resultsPage.html')
				self.send_response(200)
				self.send_header('Content-type', 'text/html')
				self.end_headers()
				self.wfile.write(f.read())
				f.close()
				return
			
			else:
				logging.debug('Incorrect PageRequest: ' + self.path)

			return

		except IOError:
			self.send_error(404,'File Not Found: %s' % self.path)
     

	def do_POST(self):
		
		logging.debug('Recieved POST request')
		global rootnode
		#try:
		ctype, pdict = cgi.parse_header(self.headers.getheader('content-type'))
		
		if ctype == 'multipart/form-data':
			logging.debug('Correct Encoding Found')
			query=cgi.parse_multipart(self.rfile, pdict)
		else:
			logging.debug('Incorrect Encoding Found')
			self.wfile.write("<HTML>POST BAD.<BR><BR>");
			return
		
		logging.debug(str(query.keys()))
		self.send_response(301)
		
		self.end_headers()
		
		temp_dir = tempfile.mkdtemp(dir = os.environ['PYTHONSCRATCH'])
		
		fasta_name = MakeFname([temp_dir, 'SeqFile.fasta'])
		ref_name = MakeFname([temp_dir, 'RefFile.gb'])
		miRNA_name = MakeFname([temp_dir, 'MiFile.txt'])
		config_name = MakeFname([temp_dir, 'config.cfg'])
		
		
		logging.debug('Writting FASTA data')
		with open(fasta_name, mode = 'w') as handle:
			if len(query['ViralSeqs_file'][0]) > 0:
				handle.write(query['ViralSeqs_file'][0])
			else:
				handle.write(query['ViralSeqs'][0])
		
		logging.debug('Checking FASTA')
		if not(CheckInputFormat(fasta_name, 'fasta')):
			logging.debug('BAD FASTA FORMAT FOUND!')
			return
		
		logging.debug('Writting GENBANK data')
		with open(ref_name, mode = 'w') as handle:
			if len(query['RefSeq_file'][0]) > 0:
				handle.write(query['RefSeq_file'][0])
			else:
				handle.write(query['RefSeq'][0])
		
		logging.debug('Checking GENBANK')
		if not(CheckInputFormat(ref_name, 'genbank')):
			logging.debug('BAD GENBANK FORMAT FOUND!')
			return
		
		
		
		logging.debug('Writting MIRNA data')
		with open(miRNA_name, mode = 'w') as handle:
			handle.write(query['InputMiRNAs'][0])
		
		this_config = ConfigParser.ConfigParser(DEFUALTCONFIG)
		this_config.set('DEFAULT', 'cache_dir', temp_dir + os.sep)
		this_config.set('DEFAULT', 'out_direc', temp_dir + os.sep)
		this_config.set('DEFAULT', 'log_name', 'log.log')
		this_config.set('DEFAULT', 'input_fasta', fasta_name)
		this_config.set('DEFAULT', 'ref_direc', temp_dir + os.sep)
		this_config.set('DEFAULT', 'scratch_dir', temp_dir + os.sep)
		this_config.set('DEFAULT', 'out_config', MakeFname([temp_dir, 'outdata.cfg']))
		this_config.set('DEFAULT', 'sub_figure', MakeFname([temp_dir, 'subtype_pie.png']))
		
		logging.debug('Writting Config data')
		with open(config_name, mode = 'w') as handle:
			this_config.write(handle)
		
		command = 'python '
		command += '"' + os.environ['MYDOCPATH'] + 'PyELM\\'
		command += 'VisualGenome.py" -v '
		command += '--config_file=' + config_name
		
		sys_call = subprocess.Popen(command)
		
		
		
		self.wfile.write('<HTML>' + 'PROCESSING' + '<BR><BR></HTML>')
		

		#except :
			#logging.debug('Error in POST request.')
			#pass

def MakeOutput(GENE_LIST):
	"""
	Creates a HTML formated output of a GENE_LIST
	@param	GENE_LIST	A SET of LLIDs
	@return				An HTML formated string
	"""
	
	output_str = string.join(map(lambda x: str(x), GENE_LIST), ',')
	
	return output_str

def CheckInputFormat(INPUT_FILE, BIO_SEQ_FORMAT):
	"""
	A utility function for making sure the input data is in the proper format.
	@param: INPUT_FILE		A file name to the user's input data
	@param: BIO_SEQ_FORMAT	The expected BioSeqFormat
	@returns:				True if the format is correct and False otherwise
	"""
	
	counter = 0
	with open(INPUT_FILE) as handle:
		if BIO_SEQ_FORMAT == 'fasta':
			for this_seq in SeqIO.parse(handle, 'fasta'):
				counter += 1
		if BIO_SEQ_FORMAT == 'genbank':
			for this_seq in SeqIO.parse(handle, 'genbank'):
				counter += 1
	
	
	if counter > 0:
		return True
	else:
		return False

def UpdateResults(TEMP_DIR):
	"""
	A thread worker which will monitor the results directory and update
	the ResultsPage.html as results are produced.
	
	@param: TEMP_DIR		The FULL_PATH to the directory where all of the 
							results will be placed.
	"""
	
	template_name = DIREC + 'ResultsPage_template.html'
	result_page = MakeFname([TEMP_DIR, 'ResultsPage.html'])
	
	num_seq_re = re.compile('XNUM_SEQSX')
	pie_chart_re = re.compile('XPIECHARTX')
	
	
	with open(template_name) as handle:
		result_string = handle.read()
	
	with open(result_page, mode = 'w') as handle:
		handle.write(result_string)
	
	out_config = ConfigParser.ConfigParser()
	with open(TEMP_DIR + os.sep + 'outdata.cfg') as handle:
		out_config.read(handle)
	
	if out_config.has_option('DEFAULT', 'wanted_seqs') & (num_seq_re != None):
		result_string = num_seq_re.sub(str(out_config.getint('DEFAULT',
									'wanted_seqs')), result_string)
		with open(result_page, mode = 'w') as handle:
			handle.write(result_string)
		num_seq_re = None
		
		
	if out_config.has_option('DEFAULT', 'sub_figure') & (pie_chart_re != None):
		result_string = pie_chart_re.sub(str(out_config.get('DEFATUL',
									'sub_figure')), result_string)
		with open(result_page, mode = 'w') as handle:
			handle.write(result_string)
			
		pie_chart_re = None
		
		
	
def main():
	
	try:
		server = HTTPServer(('', 80), MyHandler)
		print 'started httpserver...'
		server.serve_forever()
	except KeyboardInterrupt:
		print '^C received, shutting down server'
		server.socket.close()

if __name__ == '__main__':
	DIREC = os.environ['MYDOCPATH'] + 'PyELM\\'
	DEFUALTCONFIG = {'cache_dir':None,
					'out_direc':None,
					'num_threads':10,
					'verbose':False,
					'log_name':None,
					'input_fasta':None,
					'ref_direc':None,
					'wanted_subs':None,
					'elm_flag':False,
					'tf_flag':False,
					'mirna_flag':False,
					'hom_flag':False,
					'hom_file':None,
					'blast_dir':os.environ['BLASTPATH'],
					'hybrid_dir':os.environ['HYBRIDPATH'],
					'hybrid_calib': None,
					'scratch_dir': None}
	
	main()

