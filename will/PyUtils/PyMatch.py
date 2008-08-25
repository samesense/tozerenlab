from __future__ import with_statement
import re
import PyMozilla

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
	
	moz_emu = PyMozilla.MozillaEmulator(cacher = None)
	
	page_dl = moz_emu.download(search_url, search_post)
	
	post = search_re.findall(page_dl)
	seq_url = base_url + 'entrez/viewer.fcgi'
	all_seq = ''
	for this_post in post:
		post_data = 'val=%(seq)s&from=%(from)s&to=%(to)s&fmt=fasta' % \
			{'seq':this_post[0], 'from':this_post[1],'to':this_post[2]}
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
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	