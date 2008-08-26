from __future__ import with_statement
import re
import PyMozilla
from collections import defaultdict
from optparse import OptionParser
import logging

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
	
def CheckPositions(CHROM_DICT, WIDTH):
	"""
	Checks the positions defined in CHROM_DICT for TF binding positions.
	
	@param CHROM_DICT:		A dict with (key,value) as (CHR_NUM, list(POS))
	
	@returns:				A list of tuples (CHR_NUM, POS, DIST_to_nearest_TF)
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
			
			best_pos = min(tf_data, key = pos_fun)
			
			output_data.append([this_chrom, this_pos, WIDTH - best_pos[0], best_pos[5]])
			
			
	return output_data
	
if __name__ == '__main__':
	
	parser = OptionParser()
	
	parser.add_option('-f', '--filename',
						dest = 'out_file',
						type = 'string',
						action = 'store',
						help = 'Input FileName')
	parser.add_option('-v', '--verbose',
						dest = 'verbose',
						action = 'store_true',
						default = true,
						help = 'Be Verbose')
	parser.add_option('-l', '--logname',
						dest = 'log_name'
						type = 'string',
						action = 'store'
						help = 'Log FileName')
	parser.add_option('-w', '--width',
						dest = 'width'
						type = 'int',
						action = 'store'
						default = 500
						help = 'Log FileName')
	
	
	
	
	
	
	
	
	