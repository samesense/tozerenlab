import nose.tools
import PyPFAM
import threading, shelve, os








def ExampleData(NUM = 0):
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
		
	
	if NUM == 0:
		return EXAM_SEQ, EXAM_RES
	
	
	EX_2 = 'MGRITEDLIRRNAEHNDCVIFSLEELSLHQQEIERLEHIDKWCRDLKILYLQNNLIGKIENVSKLKKLEY' + \
			'LNLALNNIERIENLEGCEWLTKLDLTVNFIGELSSVKTLTHNIHLKELFLMGNPCADFDGYRQFVVVTLQ' + \
			'QLKWLDGKEIERSERIQALQNYTSVEQQIREQEKAYCLRRAKEKEEAQRKLEEENESEDKKKSSTGFDGH' + \
			'WYTDIHTACPSATENQDYPQVPETQEEQHNTKESDDIEDDLAFWNKPSLFTPESRLETLRHMEKQRKAQD' + \
			'KLSEKKKKAKPPRTLITEDGKVLNVNEAKLDFSLKDDEKHNQIILDLAVYRYMDTSLIEVDVQPTYVRVM' + \
			'VKGKPFQLALSTEVQPDRSSAKRSQTTGHLLICMPKVGEMITGGQRTPTSVKTTSTSSREQTNPRKKQIE' + \
			'RLEVDPSKHSCPDVSTIVQEKRHRPKRMESQPRDEPSEEDPDFEDNPEVPPLI'
	
	EX_2_RES = []
	
	if NUM == 1:
		return EX_2, EX_2_RES
	
	EX_3 = 'MAYRSCVVGFSSLSGCEMTPAGSPQPGTSGWGSCGLPGPGFSSRSLTSCRPAGTIPKVTVNPSLLVPLDL' + \
			'KVDPAVQQQKNQEKEEMKALNDKFASLIGKVQALEQRNQLLETRWSFLQGQGSATFDLSHHYETFQGRLQ' + \
			'EELRKVSQERGQLEANLLQVLEKVEEFRVRYEDEISKRTDLEFTFVQLKKDLDAECLRRTELETKLKGLQ' + \
			'GFLELMRTVYEQELKDLTAQVKDVSVTVGLDSRCHIDLSGIVEEVKAQYDAIAARSLEEAEAYSRSQLEE' + \
			'RAARSAEFGNSLQSSRCEIADLNVRIQKLRSQIVSVKSHCLKLEENIKVAEEQGELAFQDAKDKMAQLEN' + \
			'ALQKAKQDMARQLREYQDLMNTKLALDIEIATYHKLMEGEESRMDLPSATVVSTVKSGCRTTASKSGLTK' + \
			'TSSRKKKNRRGPVIKITEMSEKYLSQESEASE'
	
	EX_3_RES = [('Intermediate filament protein', 'Family', '82', '393', '333.2', '5.2e-97')]
	
	if NUM == 2:
		return EX_3, EX_3_RES
	
	
	
def testCheckPFAM():
	"""
	Check the BASIC PFAM
	"""
	
	for i in range(2):
		yield ChecktestCheckPFAM, i
	
	
	

def ChecktestCheckPFAM(IN):
	EXAM_SEQ, EXAM_RES = ExampleData(NUM = IN)
	
	out = PyPFAM.CheckPFAM(EXAM_SEQ)
	
	nose.tools.assert_true(out == EXAM_RES, 
			'The CheckPFAM did not produce the correct output!')
	
	
	
	
def testCheckMultiPFAM():
	"""
	Check the LIST PFAM
	"""
	
	for i in range(2):
		yield ChecktestCheckMultiPFAM, i
	
	
def ChecktestCheckMultiPFAM(IN):
	
	EXAM_SEQ, EXAM_RES = ExampleData(NUM = IN)
	
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
	
	nose.tools.assert_true(missing_out[0] == EXAM_RES[0],
		'With FORCE_LEVEL = 2 CheckPFAM should return the webserver results')
	
def testShelve():
	"""
	Test storing data in a shelve object.
	"""
	
	SHELVE_FILE = os.environ['PYTHONSCRATCH'] + 'test_shelve.slf'
	
	shelve_obj = shelve.open(SHELVE_FILE, writeback = True)
	dict_lock = threading.Lock()
	EXAM_SEQ, EXAM_RES = ExampleData()
	
	
	present_out = PyPFAM.CheckPFAM(EXAM_SEQ, DICT_OBJ = shelve_obj, 
							DICT_LOCK = dict_lock)
	
	shelve_obj.close()
	
	shelve_obj = shelve.open(SHELVE_FILE, writeback = True)
	
	this_out = shelve_obj[EXAM_SEQ]
	
	nose.tools.assert_true(this_out == EXAM_RES, 
			'The shelve did not keep the correct output')
	
	shelve_obj.close()
	
	
	os.remove(SHELVE_FILE)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	