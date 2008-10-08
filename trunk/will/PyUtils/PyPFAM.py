from __future__ import with_statement
import re, threading, time, logging
import PyMozilla



def CheckMultiPFAM(INPUT_SEQ_LIST, NUM_THREADS = 10, DICT_OBJ = None, 
					DICT_LOCK = None, SHELVE_FILE = None, FORCE_LEVEL = 1):
	"""
	Uses multithreading to process multiple requests at a time.  Overloading 
	the server by using NUM_THREADS > 10 may not increase the speed due to 
	limiting factors on the SERVER side architechture.
	
	@param:	INPUT_SEQ_LIST		A LIST of protein sequence to be submitted 
								to the PFAM webserver.  
								No checking is performed!!!
	
	@param:	DICT_OBJ = None	A python shelve object which is key'ed by the 
								INPUT_SEQ and contains the result of a previous
								CheckPFAM run.
	
	@param:	DICT_LOCK = None	A threading Lock object which is shared between
								multiple threads and prevents simultatious 
								read/writes that will corrupt the cache.
	
	@param:	NUM_THREADS = 10	The number of concurrent threads to pound the 
								server with.
	
	
	This will return a list of lists of tuples as described by CheckPFAM
	"""
	
	def PFAMWorker(INDEX, SEQ):
		output = CheckPFAM(SEQ, DICT_OBJ = DICT_OBJ, DICT_LOCK = DICT_LOCK,
							FORCE_LEVEL = FORCE_LEVEL)
		with list_lock:
			output_list[INDEX] = output
		worker_seph.release()
	
	
	if DICT_OBJ == None:
		DICT_OBJ = {}
		DICT_LOCK = threading.Lock()
	
	
	worker_seph = threading.Semaphore(NUM_THREADS)
	list_lock = threading.Lock()
	
	#initialize it to numbers just so we have a place to put answers when we get them
	output_list = range(len(INPUT_SEQ_LIST))
	
	all_threads = []
	
	for this_ind in xrange(len(INPUT_SEQ_LIST)):
		logging.debug('Waiting for free Thread')
		time.sleep(5)
		worker_seph.acquire()
		
		this_thread = threading.Thread(target = PFAMWorker, 
							args = (this_ind, 
									INPUT_SEQ_LIST[this_ind]))
		this_thread.start()
		all_threads.append(this_thread)
		
	for this_thread in all_threads:
		this_thread.join()
	
	#after all threads have been .join()ed then all spots have been filled
	
	
	return output_list
	
	

def CheckPFAM(INPUT_SEQ, DICT_OBJ = None, DICT_LOCK = None, 
				FORCE_LEVEL = 1):
	"""
	Checks the INPUT_SEQ against the PFAM database webserver.
	returns a list of tuples:
	
	@param:	INPUT_SEQ			A protein sequence to be submitted to the PFAM
								webserver.  No checking is performed!!!
	
	@param:	DICT_OBJ = None	A python shelve object which is key'ed by the 
								INPUT_SEQ and contains the result of a previous
								CheckPFAM run.
	
	@param:	DICT_LOCK = None	A threading Lock object which is shared between
								multiple threads and prevents simultatious 
								read/writes that will corrupt the cache.
	
	@param:	FORCE_LEVEL = (int)	A choice variable to control whether to fetch 
								results or not.
	
		FORCE_LEVEL = 0			Return a cached result if present and return 
								None if the result is not cached ... useful
								for quick checks.
		
		FORCE_LEVEL = 1			Return a cached result if present and otherwise
								perform a webquerry and cache the results. 
								(DEAFULT)
		
		FORCE_LEVEL = 2			ALWAYS perform a webquery and then cache the 
								results.  Useful for when the database updates.
	
	@returns:  A list of tuples:
	
	[(NAME_1, TYPE_1, START_POS_1, END_POS_1, BIT_SCORE_1, E_VAL_1),
	 (NAME_2, TYPE_2, START_POS_2, END_POS_2, BIT_SCORE_2, E_VAL_1),
	 ...]
	 
	 ALL VALUES ARE STRINGS!!!!!!!!!!!!!!
	"""
	
	#if a shelve object is provided then we MUST have a Lock for it!
	if (DICT_OBJ != None) & (DICT_LOCK == None):
		raise TypeError, 'If DICT_FILE is provided then DICT_LOCK ' \
						+ 'cannot be None'
	
	#check to see whether we've already found this sequence
	if DICT_OBJ != None:
		with DICT_LOCK:
			cached = DICT_OBJ.has_key(INPUT_SEQ)
		if cached:
			if FORCE_LEVEL <= 1:
				logging.debug('Returning cached object')
				return DICT_OBJ[INPUT_SEQ]
		else:
			if FORCE_LEVEL == 0:
				logging.debug('Returning None')
				return None
	
	BASE_URL = 'http://pfam.sanger.ac.uk/search/sequence'
	
	moz_emu = PyMozilla.MozillaEmulator(cacher=None)
	
	logging.debug('Submitting SEQ to webserver')
	splash_page = moz_emu.post_multipart(BASE_URL, [('seq', INPUT_SEQ)], (), 
					trycount = 2)
	
	time_out, stat_url, job_id, result_url = ExtractCheckingInfo(splash_page)
	
	while time_out != 0:
		logging.debug('Sleeping for %(time)d' % {'time':time_out})
		time.sleep(time_out)
		time_out = CheckStatus(moz_emu, stat_url, job_id)
	
	
	final_output = ExtractINFO(moz_emu, result_url, job_id)
	
	if DICT_OBJ != None:
		logging.debug('Caching results for later use')
		with DICT_LOCK:
			DICT_OBJ[INPUT_SEQ] = final_output
	
	return final_output
	
	
	
def ExtractCheckingInfo(SPLASH_PAGE):
	"""
	Uses regular expressions to extract the CheckURL, JOB_ID, 
	ESTIMATED_TIME, and RESULT_URL.
	
	@param: SPLASH_PAGE			The page returned from submitting the sequence
								to the PFAM webserver.
								
	@returns:
	EST_TIME					The estimated time (as an int) for completion
	STAT_URL					The URL to check for the status of the job
	JOB_ID						The jobId of this submission
	RESULT_URL					The URL to check when the results are finished.
	"""
	
	reg_str = '"estimatedTime":(\d*),'
	reg_str += '"checkURI":"(.*?)",'
	reg_str += '.*?"jobId":"(.*?)",'
	reg_str += '.*?,"doneURI":"(.*?)"'
	
	output = re.findall(reg_str, SPLASH_PAGE)
	
	if len(output) != 1:
		bad_str = 'The SPLASH_PAGE returned to many URLs to check' +str(output)
		raise IndexError, bad_str
	
	if len(output[0]) != 4:
		bad_str = 'The SPLASH_PAGE did not have all of the information '
		bad_str += 'required' + str(output)
		raise IndexError, bad_str
	
	est_time = int(output[0][0])
	stat_url = output[0][1]
	job_id = output[0][2]
	result_url = output[0][3]
	
	return est_time, stat_url, job_id, result_url
	
	
def CheckStatus(MOZ_EMU, CHECK_URL, JOB_ID):
	"""
	Checks the CHECK_URL and returns a wait time.  Will return 0 if the 
	status is "DONE"
	
	@param:	MOZ_EMU			A PyMozilla Emulator object.
	
	@param:	CHECK_URL		The URL to check the status at.
	
	@param:	JOB_ID			The jobId of the submitted sequence.
	
	@returns:
		int(TIME)			The time to wait before checking again.
							Returns 0 if the job has finished!
	"""
	
	logging.debug("Checking for DONE'ness")
	check_out = MOZ_EMU.post_multipart(CHECK_URL, [('jobId', JOB_ID)],(), 
										trycount = 2)
	output = re.findall('"status":"(.*?)"', check_out)
	
	
	if output[0] == 'PEND':
		logging.debug('Job Pending')
		return 10
	elif output[0] == 'RUN':
		logging.debug('Job Running')
		return 10
	elif output[0] == 'DONE':
		logging.debug('Webserver finished processing')
		return 0
	else:
		return 10
	
	
def ExtractINFO(MOZ_EMU, RESULT_URL, JOB_ID):
	"""
	Extract the PFAM information from the website.
	
	@param:	MOZ_EMU			A PyMozilla Emulator object.
	
	@param:	RESULT_URL		A URL which will have the result information.
	
	@param:	JOB_ID			The jobId of the submitted sequence.
	
	@returns:
	A list of tuples:
	
	[(NAME_1, TYPE_1, START_POS_1, END_POS_1, BIT_SCORE_1, E_VAL_1),
	 (NAME_2, TYPE_2, START_POS_2, END_POS_2, BIT_SCORE_2, E_VAL_1),
	 ...]
	 
	 ALL VALUES ARE STRINGS!!!!!!!!!!!!!!
	
	"""
	
	reg_str = '<tr class=".*?Significant">.*?' #only find significant classes
	reg_str += '<td class="desc">(.*?)</td>.*?' #catch the PFAM name
	reg_str += '<td>(.*?)</td>.*?' #catch the type
	reg_str += '<td>(\d*)</td>.*?' #catch the start_pos
	reg_str += '<td>(\d*)</td>.*?' #catch the end_pos
	reg_str += '<td>.*?</td>.*?' #discard HMM start_pos
	reg_str += '<td>.*?</td>.*?' #discard HMM end_pos
	reg_str += '<td>(.*?)</td>.*?' #catch bit-score
	reg_str += '<td>(.*?)</td>' #catch e-value
	
	re_obj = re.compile(reg_str, re.S)
	
	logging.debug('Getting Results Page')
	res_page = MOZ_EMU.post_multipart(RESULT_URL, [('jobId', JOB_ID)], (),
										trycount = 2)
	
	output = re_obj.findall(res_page)
	
	return output
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	


	
	