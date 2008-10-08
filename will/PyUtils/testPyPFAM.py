import nose.tools
import PyPFAM
import threading








def ExampleData():
	EXAM_SEQ = 'MAGAASPCANGCGPSAPSDAEVVHLCRSLEVGTVMTLFYSKKSQRPERKTFQVKLETRQI' + \
	'TWSRGADKIEGAIDIREIKEIRPGKTSRDFDRYQEDPAFRPDQSHCFVILYGMEFRLKTL' + \
	'SLQATSEDEVNMWIRGLTWLMEDTLQAATPLQIERWLRKQFYSVDRNREDRISAKDLKNM' + \
	'LSQVNYRVPNMRFLRERLTDLEQRTSDITYGQFAQLYRSLMYSAQKTMDLPFLEASALRA' + \
	'GERPELCRVSLPEFQQFLLEYQGELWAVDRLQVQEFMLSFLRDPLREIEEPYFFLDEFVT' + \
	'FLFSKENSIWNSQLDEVCPDTMNNPLSHYWISSSHNTYLTGDQFSSESSLEAYARCLRMG' + \
	'CRCIELDCWDGPDGMPVIYHGHTLTTKIKFSDVLHTIKEHAFVASEYPVILSIEDHCSIA' + \
	'QQRNMAQYFKKVLGDTLLTKPVDIAADGLPSPNQLKRKILIKHKKLAEGSAYEEVPTSVM' + \
	'YSENDISNSIKNGILYLEDPVNHEWYPHYFVLTSSKIYYSEETSSDQGNEDEEEPKEASG' + \
	'STELHSNEKWFHGKLGAGRDGRHIAERLLTEYCIETGAPDGSFLVRESETFVGDYTLSFW' + \
	'RNGKVQHCRIHSRQDAGTPKFFLTDNLVFDSLYDLITHYQQVPLRCNEFEMRLSEPVPQT' + \
	'NAHESKEWYHASLTRAQAEHMLMRVPRDGAFLVRKRNEPNSYAISFRAEGKIKHCRVQQE' + \
	'GQTVMLGNSEFDSLVDLISYYEKHPLYRKMKLRYPINEEALEKIGTAEPDYGALYEGRNP' + \
	'GFYVEANPMPTFKCAVKALFDYKAQREDELTFTKSAIIQNVEKQEGGWWRGDYGGKKQLW' + \
	'FPSNYVEEMVSPAALEPEREHLDENSPLGDLLRGVLDVPACQIAVRPEGKNNRLFVFSIS' + \
	'MASVAHWSLDVAADSQEELQDWVKKIREVAQTADARLTEGKMMERRKKIALELSELVVYC' + \
	'RPVPFDEEKIGTERACYRDMSSFPETKAEKYVNKAKGKKFLQYNRLQLSRIYPKGQRLDS' + \
	'SNYDPLPMWICGSQLVALNFQTPDKPMQMNQALFLAGGHCGYVLQPSVMRDEAFDPFDKS' + \
	'SLRGLEPCAICIEVLGARHLPKNGRGIVCPFVEIEVAGAEYDSIKQKTEFVVDNGLNPVW' + \
	'PAKPFHFQISNPEFAFLRFVVYEEDMFSDQNFLAQATFPVKGLKTGYRAVPLKNNYSEGL' + \
	'ELASLLVKIDVFPAKQENGDLSPFGGASLRERSCDASGPLFHGRAREGSFEARYQQPFED' + \
	'FRISQEHLADHFDGRDRRTPRRTRVNGDNRL'
	
	EXAM_RES = [('PH domain', 'Domain', '33', '142', '40.8', '6.1e-11'), 
		('Phosphatidylinositol-specific phospholipase C, X domain', 'Family', 
		'322', '465', '328.4', '1.4e-95'), 
		('SH2 domain', 'Domain', '550', '639', '3.5e-32', 'ls'), 
		('SH2 domain', 'Domain', '668', '741', '1.7e-29', 'ls'), 
		('SH3 domain', 'Domain', '794', '849', '85.6', '1.7e-22'),
		('Phosphatidylinositol-specific phospholipase C, Y domain', 'Family', 
		'952', '1070', '177.2', '4.6e-50'), 
		('C2 domain', 'Domain', '1090', '1177', '73.7', '3.2e-21')]
	
	return EXAM_SEQ, EXAM_RES
	
def testCheckPFAM():
	"""
	Check the BASIC PFAM
	"""
	
	EXAM_SEQ, EXAM_RES = ExampleData()
	
	out = PyPFAM.CheckPFAM(EXAM_SEQ)
	
	nose.tools.assert_true(out == EXAM_RES, 
			'The CheckPFAM did not produce the correct output!')
	
	
	
def testCheckMultiPFAM():
	"""
	Check the LIST PFAM
	"""
	
	EXAM_SEQ, EXAM_RES = ExampleData()
	
	test_list = [EXAM_SEQ, EXAM_SEQ, EXAM_SEQ, EXAM_SEQ]
	
	out_list = PyPFAM.CheckMultiPFAM(test_list)
	
	for this_check in out_list:
		nose.tools.assert_true(this_check == EXAM_RES, 
			'The CheckPFAM did not produce the correct output!')
			
def testFLEVEL0():
	"""
	Test the FORCE_LEVEL = 0
	"""
	
	EXAM_SEQ, EXAM_RES = ExampleData()
	
	test_dict = {EXAM_SEQ:EXAM_RES}
	#needed for the 
	junk_lock = threading.Lock()
	
	present_out = PyPFAM.CheckPFAM(EXAM_SEQ, DICT_OBJ = test_dict, 
							DICT_LOCK = junk_lock, 
							FORCE_LEVEL = 0)
	
	nose.tools.assert_true(present_out == EXAM_RES, 
			'CheckPFAM did not return the correct result')
			
	missing_out = PyPFAM.CheckPFAM(EXAM_SEQ[0:-1], DICT_OBJ = test_dict, 
							DICT_LOCK = junk_lock, 
							FORCE_LEVEL = 0)
	
	nose.tools.assert_true(missing_out == None,
	'With FORCE_LEVEL = 0 CheckPFAM should return None for missing elements')
	
def testFLEVEL1():
	"""
	Test the FORCE_LEVEL = 1
	"""
	
	EXAM_SEQ, EXAM_RES = ExampleData()
	
	test_dict = {EXAM_SEQ:EXAM_RES, EXAM_SEQ[1:]: 'JUNK_DATA'}
	#needed for the 
	junk_lock = threading.Lock()
	
	present_out = PyPFAM.CheckPFAM(EXAM_SEQ, DICT_OBJ = test_dict, 
							DICT_LOCK = junk_lock, 
							FORCE_LEVEL = 1)
	
	nose.tools.assert_true(present_out == EXAM_RES, 
			'CheckPFAM did not return the correct result')
			
	missing_out = PyPFAM.CheckPFAM(EXAM_SEQ[0:-1], DICT_OBJ = test_dict, 
							DICT_LOCK = junk_lock, 
							FORCE_LEVEL = 1)
	
	nose.tools.assert_true(missing_out == EXAM_RES,
		'With FORCE_LEVEL = 1 CheckPFAM should return the fetched results')
	
	missing_out = PyPFAM.CheckPFAM(EXAM_SEQ[1:], DICT_OBJ = test_dict, 
							DICT_LOCK = junk_lock, 
							FORCE_LEVEL = 1)
	
	nose.tools.assert_true(missing_out == 'JUNK_DATA',
		'With FORCE_LEVEL = 1 CheckPFAM should return the cached results')
	
def testFLEVEL2():
	"""
	Test the FORCE_LEVEL = 2
	"""
	
	EXAM_SEQ, EXAM_RES = ExampleData()
	
	test_dict = {EXAM_SEQ:EXAM_RES, EXAM_SEQ[1:]: 'JUNK_DATA'}
	#needed for the 
	junk_lock = threading.Lock()
	
	present_out = PyPFAM.CheckPFAM(EXAM_SEQ, DICT_OBJ = test_dict, 
							DICT_LOCK = junk_lock, 
							FORCE_LEVEL = 2)
	
	nose.tools.assert_true(present_out == EXAM_RES, 
			'CheckPFAM did not return the correct result')
			
	missing_out = PyPFAM.CheckPFAM(EXAM_SEQ[0:-1], DICT_OBJ = test_dict, 
							DICT_LOCK = junk_lock, 
							FORCE_LEVEL = 2)
	
	nose.tools.assert_true(missing_out == EXAM_RES,
		'With FORCE_LEVEL = 2 CheckPFAM should return the fetched results')
	
	missing_out = PyPFAM.CheckPFAM(EXAM_SEQ[0:-2], DICT_OBJ = test_dict, 
							DICT_LOCK = junk_lock, 
							FORCE_LEVEL = 2)
	s
	nose.tools.assert_true(missing_out[0] == EXAM_RES[0],
		'With FORCE_LEVEL = 2 CheckPFAM should return the webserver results')
	
	
	
	
	
	
	
	
	