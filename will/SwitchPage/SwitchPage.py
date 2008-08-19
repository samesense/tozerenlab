from __future__ import with_statement
import string,cgi,time,os,logging,re, PyMozilla

from BaseHTTPServer import BaseHTTPRequestHandler, HTTPServer

logging.basicConfig(level=logging.DEBUG,
					format='%(asctime)s %(levelname)s %(message)s')
					
class MyHandler(BaseHTTPRequestHandler):

	def do_GET(self):
		logging.debug('Recieved GET request')
		try:
			if self.path.endswith('SwitchPage.htm'):
				logging.debug('Correct PageName')
				f = open(DIREC + 'SwitchPage.htm') #self.path has /test.html
#note that this potentially makes every file on your computer readable by the internet

				self.send_response(200)
				self.send_header('Content-type',	'text/html')
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
			
		self.send_response(301)
		
		self.end_headers()
		
		norm_list = CheckGeneList(query['GeneList'])
		
		these_switches = norm_list.intersection(SWITCH_GENES)
		
		output_str = MakeOutput(these_switches)
		
		self.wfile.write('<HTML>' + output_str + '<BR><BR></HTML>')
		

		#except :
			#logging.debug('Error in POST request.')
			#pass

def CheckGeneList(GENE_LIST):
	"""
	Checks the provided gene list and returns the gene that are in the 
	database and thier LLIDs.
	@param:	GENE_LIST	A list of input genes.
	@return:			A set of LLIDs.
	"""
	
	output_list = set()
	
	for this_gene in GENE_LIST[0].split('\r'):
		stripped_gene = this_gene.strip()
		
		if CONVERSION_DICT.has_key(stripped_gene):
			output_list.add(CONVERSION_DICT[stripped_gene])
			out_str = 'Found gene: '
			out_str += stripped_gene + ' as '
			out_str += str(CONVERSION_DICT[stripped_gene])
			logging.debug(out_str)
		else:
			logging.debug('Missing gene: ' + stripped_gene)
	return output_list
	
	
def RemakeSwitchPage():
	"""
	Creates the SwitchPage.htm based on the phenotype.txt file
	"""
	
	TAG = '\t%%%%TemplatePhen%%%%'
	
	with open(DIREC + 'SwitchPage_TEMPLATE.htm') as handle:
		page_data = handle.read()
	
	ind = page_data.find(TAG)
	
	out_page_data = page_data[0:ind]
	
	with open(DIREC + 'phenotypes.txt') as handle:
		for this_phen in handle:
			out_page_data += '\t<option>' 
			out_page_data += this_phen[0:-1]
			out_page_data += '</option>\n' 
	
	out_page_data += page_data[ind+len(TAG):]
	
	with open(DIREC + 'SwitchPage.htm', mode = 'w') as handle:
		handle.write(out_page_data)
	
def MakeTransDict(FILE_NAME):
	"""
	Loads a tab-delimited file and makes a dictionary which maps the IDs back 
	to LLIDs.
	@param:		FILE_NAME		The path to the gene_IDS.txt
	@return:	map_to_LLID		A dictionary which maps Entrez IDs, 
								Gene Symbols and Affy IDS to LLIDs
	"""
	
	map_to_LLID = {}
	
	with open(FILE_NAME) as handle:
		for this_line in handle:
			parts = this_line.split('\t')
			if len(parts[0]) != 0:
				this_llid = int(parts[0])
				map_to_LLID[parts[0].strip()] = this_llid
			else:
				continue
			if len(parts[1]) != 0:
				map_to_LLID[parts[1].strip()] = this_llid
			if len(parts[2]) != 0:
				map_to_LLID[parts[2].strip()] = this_llid
	
	return map_to_LLID

def LoadSwitchGenes(FILE_NAME):
	"""
	Loads the switch_genes.txt file into a SET of integers.
	"""
	
	switch_genes = set()
	
	with open(FILE_NAME) as handle:
		for this_line in handle:
			if len(this_line) != 0:
				switch_genes.add(int(this_line))
	
	return switch_genes
	
def MakeOutput(GENE_LIST):
	"""
	Creates a HTML formated output of a GENE_LIST
	@param	GENE_LIST	A SET of LLIDs
	@return				An HTML formated string
	"""
	
	output_str = string.join(map(lambda x: str(x), GENE_LIST), ',')
	
	return output_str
	
	
def main():
	RemakeSwitchPage()

	try:
		server = HTTPServer(('', 80), MyHandler)
		print 'started httpserver...'
		server.serve_forever()
	except KeyboardInterrupt:
		print '^C received, shutting down server'
		server.socket.close()

if __name__ == '__main__':
	DIREC = os.environ['MYDOCPATH'] + 'SwitchPage\\'
	
	CONVERSION_DICT = MakeTransDict(DIREC + 'gene_IDS.txt')
	SWITCH_GENES = LoadSwitchGenes(DIREC + 'switch_genes.txt')
	
	
	
	main()

