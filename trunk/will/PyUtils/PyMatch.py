from __future__ import with_statement
import re
import PyMozilla
from collections import defaultdict
from optparse import OptionParser, OptionGroup
import logging
import string

def TFChecker(INPUT_SEQ):
	"""
	Uses a MATCH login to find the transcription factor binding sites in 
	INPUT_SEQ.
	
	Returns a list of lists:
		[[pos1, strand1, core_sim1, total_sim1, seq1, name1, mat1],
		 [[pos2, strand2, core_sim2, total_sim2, seq2, name2, mat2]]
	
	"""
	moz_emu = PyMozilla.MozillaEmulator(cacher=None)
	info_getter = re.compile('<A HREF.*?</a>(.*?)<a href.*?</a>',re.S)
	mat_getter = re.compile('<A HREF.*?>(.*?)</a>')
	name_getter = re.compile('<a href.*?>(.*?)</a>')
	sep = re.compile(' *')
	base_url = 'http://www.gene-regulation.com'
	user_name = 'walldo2'
	password = 'bfg20k'
	login_post = 'user=' + user_name
	login_post += '&password=' + password
	login_post += '&request_uri=%2Findex.html&Log+In=Log+In'
	
	#login and get the cookie
	logging.info('Logging into MATCH')
	moz_emu.download(base_url + '/login', login_post)
	
	match_base = '/cgi-bin/pub/programs/match/bin/match.cgi'
	post_data = 'Status=First'
	post_data += '&searchName=default'
	post_data += '&usr_seq=default.seq'
	post_data += '&seqStat=DEL'
	post_data += '&sequenceName=default'
	post_data += '&&theSequence=' + INPUT_SEQ
	post_data += '&SearchMode=ALL'
	post_data += '&group=vertebrates'
	post_data += '&quality=high'
	post_data += '&allMode=FP'
	post_data += '&cut-off=0.7'
	post_data += '&core=0.75'
	post_data += '&ourProfile=muscle_specific.prf'
	
	logging.info('Getting Match Data')
	output = moz_emu.download(base_url + match_base,post_data)
	
	val = map(string.strip,info_getter.findall(output))
	info = map(sep.split,val)
	tf_names = name_getter.findall(output)
	mat_names = mat_getter.findall(output)
	
	if len(val) == 0:
		return []
	
	converter_fun = lambda x:[int(x[0]), x[1], float(x[2]), float(x[3]), x[4]]
	
	final_output = []
	for i in range(len(info)):
		this_out = converter_fun(info[i])
		this_out.append(tf_names[i])
		this_out.append(mat_names[i])
		final_output.append(this_out)
	return final_output
	
def GetSequence(CHR_NUM, BEGIN, END):
	"""
	Retrieves the human chromosomal sequence from NCBI.
	
	@param CHR_NUM:		The chromosomal number
	@param BEGIN:		The beggining of the desired sequence
	@param END:			The end of the desired sequence region
	
	@returns:			The sequence between BEGIN and END
	"""
	
	search_re = re.compile('.*?entrez/viewer.cgi\?val=(.*?)\&from=(\d*)\&to=(\d*).*?" onCLick.*?>Save to Disk<')
	
	fasta_re = re.compile("<div class='recordbody'>(.*?)</div>", flags = 16)
	
	base_url = 'http://www.ncbi.nlm.nih.gov/'
	
	search_url = base_url + 'projects/mapview/seq_reg.cgi'
	
	search_post = 'TAXID=9606'
	search_post += '&chr=%(chrom_num)s' % {'chrom_num':str(CHR_NUM)}
	search_post += '&strand=1'
	search_post += '&from=%(begin)s' % {'begin':str(BEGIN)}
	search_post += '&to=%(end)s' %{'end':str(END)}
	search_post += '&&fmt=fasta'
	
	logging.info('Downloading Sequence:\t' + CHR_NUM + ':' + str(BEGIN))
	moz_emu = PyMozilla.MozillaEmulator(cacher = None)
	
	page_dl = moz_emu.download(search_url, search_post)
	
	post = search_re.findall(page_dl)
	seq_url = base_url + 'entrez/viewer.fcgi'
	all_seq = ''
	for this_post in post:
		post_data = 'val=%(seq)s&from=%(from)s&to=%(to)s&fmt=fasta' % \
			{'seq':this_post[0], 'from':this_post[1],'to':this_post[2]}
		
		logging.info('Retrieving Sequence:\t' + CHR_NUM + ':' + str(BEGIN))
		this_dl = moz_emu.download(seq_url, post_data)
		this_fasta = fasta_re.findall(this_dl)
		parts = this_fasta[0].split('\n')
		for this_line in parts[1:]:
			all_seq += this_line.upper()
	
	if len(all_seq) - 1 != END - BEGIN:
		num_missing = abs((END - BEGIN) - len(all_seq))
		str_dict = {'m': num_missing, 'c': CHR_NUM, 'b':BEGIN, 'e':END}
		bad_str = 'Could not retrieve %(m)d bp' % str_dict
		bad_str += ' for region %(c)s [%(b)d:%(e)d]' % str_dict
		logging.warning(bad_str)
	
	return all_seq
	
def ParseInput(FILE_NAME):
	"""
	Parses the input file to order things by chromosome and determine which
	regions to retrieve.
	
	File FMT
	
	CHR_NUM \t BASE_POS
	"""
	
	chrom_dict = defaultdict(set)
	
	
	with open(FILE_NAME) as handle:
		for this_line in handle:
			parts = this_line.split('\t')
			chrom_dict[parts[0]].add(int(parts[1]))
	
	return chrom_dict
	
def CheckPositions(CHROM_DICT, WIDTH, TOL):
	"""
	Checks the positions defined in CHROM_DICT for TF binding positions.
	
	@param CHROM_DICT:		A dict with (key,value) as (CHR_NUM, list(POS))
	
	@returns:				A list of lists [CHR_NUM, POS, DIST_to_nearest_TF]
	"""
	
	
	pos_fun = lambda x: abs(x[0] - WIDTH)
	output_data = []
	
	for this_chrom in CHROM_DICT:
		for this_pos in CHROM_DICT[this_chrom]:
			start_pos = this_pos - WIDTH
			end_post = this_pos + WIDTH
			seq_data = GetSequence(this_chrom, start_pos, end_post)
			
			tf_data = TFChecker(seq_data)
			if len(tf_data) == 0:
				continue
			
			
			tf_data.sort(key = pos_fun)
			for this_tf in tf_data:
				if pos_fun(this_tf) > TOL:
					break
				else:
					str_dict = {'c':this_chrom, 'p':this_pos, 
								'o':WIDTH - this_tf[0], 'n':this_tf[5]}
					out_str = 'FOUND TF! %(c)s:%(p)d:%(o)d:%(n)s' % str_dict
					logging.warning(out_str)
					output_data.append([this_chrom, this_pos, WIDTH - this_tf[0], this_tf[5]])
			
			
	return output_data
	
if __name__ == '__main__':
	
	parser = OptionParser(usage = 'usage: %prog [options] input_file')
	
	parser.add_option('-v', '--verbose',
						dest = 'verbose',
						action = 'store_true',
						default = False,
						help = 'Be Verbose')
	parser.add_option('-l', '--logname',
						dest = 'log_name',
						type = 'string',
						action = 'store',
						help = 'Log FileName')
	parser.add_option('-o', '--outputfile',
						dest = 'out_file',
						type = 'string',
						action = 'store',
						help = 'Output FileName')
	
	group = OptionGroup(parser, 'Sequence Options')
	
	group.add_option('-w', '--width',
						dest = 'width',
						type = 'int',
						action = 'store',
						default = 500,
						help = 'The sequence width to retrieve from NCBI')
	group.add_option('-t', '--tolerance',
						dest = 'tol',
						type = 'int',
						action = 'store',
						default = 5,
						help = 'The width to accept a TF')
	parser.add_option_group(group)
	
	(options, args) = parser.parse_args()
	
	file_fmt = '%(levelname)s \t %(threadName)s \t %(funcName)s \t %(asctime)s \t %(message)s'
	console_fmt = '%(asctime)-7s %(message)s'
	
	
	
	
	logging.basicConfig(level = logging.INFO)
	
	console = logging.StreamHandler()
	if options.verbose:
		console.setLevel(logging.INFO)
	else:
		console.setLevel(logging.WARNING)
	
	formatter = logging.Formatter(console_fmt)
	console.setFormatter(formatter)
	logging.getLogger('').addHandler(console)
	
	if options.log_name != None:
		file_log = logging.StreamHandler()
		formatter = logging.Formatter(file_fmt)
		file_log.setLevel(logging.INFO)
		file_log.setFormatter(formatter)
		logging.getLogger('').addHandler(file_log)
	
	input_data = ParseInput(args[0])
	
	output_data = CheckPositions(input_data, options.width, options.tol)
	str_conv = lambda x: str(x)
	with open(options.out_file, mode = 'w') as handle:
		for this_out in output_data:
			out_list = map(str_conv, this_out)
			out_str = string.join(out_list, '\t')
			handle.write(out_str + '\n')
	
	
	
	
	