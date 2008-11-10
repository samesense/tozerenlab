from __future__ import with_statement
import itertools as IT
import HIVGenoTypeChecker as GenoTyping
import subprocess
import tempfile
import threading
import time
import types
import string
import os
import re
import types
import Queue
import PyBLAST
import string
import numpy
import bisect
from Bio import SeqIO
from Bio import SeqFeature as SeqFeatClass
#a hack but the only way I can get everything to work right
from Bio.SeqFeature import SeqFeature
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein, generic_nucleotide
from Bio.Blast import NCBIXML, NCBIStandalone
from GenomeDiagram import *
from reportlab.lib import colors
from AnnotUtils import *
from collections import defaultdict

import logging



class ViralSeq():
	def __init__(self, TEST_SEQ, SEQ_NAME):
		#if SEQ_NAME is None then assume that TEST_SEQ is a SeqRecord
		if SEQ_NAME == None:
			self.my_sequence = TEST_SEQ.seq.tostring().upper()
			self.seq_name = TEST_SEQ.id
		else:
			self.my_sequence = TEST_SEQ
			self.seq_name = SEQ_NAME
		self.annotation = None
		self.feature_annot = []
		self.feature_annot_type = {}
		self.this_genome = None
		

	def __hash__(self):
		return hash(self.my_sequence + self.seq_name)

	def __getattr__(self, NAME):
		if NAME == 'tested_subtype':
			self.DetSubtype()
			return self.tested_subtype
		if NAME == 'id':
			return self.GenomeToBioPython().id
		if NAME == 'description':
			return self.GenomeToBioPython().description
		if NAME == 'seq':
			return self.GenomeToBioPython().seq
		if NAME == 'my_sequence':
			self.GetSequence()
			return self.my_sequence
		else:
			raise AttributeError, NAME
	
	def GenomeToBioPython(self):
		"""
		GenomeToBioPython
				Returns a generic SeqRecord in the biopython format for
		the whole genome sequence.
		"""
		return SeqRecord(Seq(self.my_sequence, generic_nucleotide),
						 id=self.seq_name)

	def AnnotToBioPython(self):
		"""
		GenomeToBioPython
				Returns a list of generic SeqRecords in the biopython
		format for each protein in self.annotation.
		"""
		final_list = []
		for this_gene in self.annotation:
			this_record = SeqRecord(Seq(self.annotation[this_gene].aa_seq,
										generic_protein),
										id = self.seq_name + ':' + this_gene)
			final_list.append(this_record)
		return final_list
	
	def SetOffset(self, REF_SEQ):
		"""
		A utility function which provides backcompatibilty with the SetOffset
		of PatSeq.  This function does nothing for RefSeqs and PyVirus's
		"""
		pass
	
	
	
	def GetSeqFeatures(self, PROT_REL = False):
		"""
		GetSeqFeatures
			Returns a list of SeqFeatures describing each protein in 
			self.annotation suitalbe for using with GenomeDiagram
		"""
		feat_list = []
		start = 0
		for this_gene in self.annotation:
			self.annotation[this_gene].GeneAnnot(None, start, 
								start + len(self.annotation[this_gene].aa_seq))
			start += len(self.annotation[this_gene].aa_seq)
			if PROT_REL:
				feat_list.append(self.annotation[this_gene].GetSeqFeature_PROT())
			else:
				feat_list.append(self.annotation[this_gene].GetSeqFeature())
		
		return feat_list

	def DetSubtype(self):
		"""
		DetSubtype
				Determines the subtype of the sequence using the
		HIVSubtype module.
		"""
		logging.debug('Determining Subtype: ' + self.seq_name)
		self.tested_subtype = GenoTyping.GetSimple(self.my_sequence)
		self.tested_subtype = self.tested_subtype.split('|')[1]
		logging.debug( self.seq_name + ': ' + self.tested_subtype)
		
	def TranslateAll(self, REFBASE, WANTED_REF = None):
		"""
		TranslateAll
				Uses a BLASTx query to determine the which ORFs
		correspond to HIV-1 proteins.  It then stores the data in a
		dictionary keyed by "genename".
		"""
		blast_record = REFBASE.BLASTx(self.GenomeToBioPython())
		gene_dict = {}
		reg = re.compile('.*?\|(.*?):(\w{3}).*')
		aa_start = 0
		for this_check in blast_record.alignments:
			label = reg.match(this_check.title).groups()
			this_gene = str(label[1])
			this_ref = str(label[0])
			if (WANTED_REF != None) & (this_ref != WANTED_REF):
				#make sure we translate "in context" of a reference
				continue
			if this_check.hsps[0].score < 900:
				continue
			if not(gene_dict.has_key(this_gene)):
				aa_seq = str(this_check.hsps[0].query)
				start_pos = this_check.hsps[0].query_start
				end_pos = start_pos + this_check.length*3
				nt_seq = self.my_sequence[start_pos:end_pos]
				gene = Gene(this_gene, this_gene.upper(),
							nt_seq, aa_seq, start_pos, end_pos)
				gene.GeneAnnot(None, aa_start, aa_start + len(aa_seq))
				
				aa_start += len(aa_seq)
				gene_dict[this_gene] = gene
			
			#end indent back
		self.annotation = gene_dict
		
		self.SetOffset(REFBASE.GetRefSeq(WANTED_REF))
		
		
	def CheckSeq(self, INPUT_SEQ):
		"""
		CheckSeq
			Returns True if the sequence is present within the sequence and False
			otherwise.
		"""
		return self.my_sequence.find(INPUT_SEQ) != -1
	
	def WriteGenesToDiagram(self, INPUT_TRACK, PROT = False):
		"""
		WriteGenomeDiagram
			Use the GenomeDiagram module to write a publication quality figure
			of the genes and features in self.annotation
		"""
		
		gene_feat_list = self.GetSeqFeatures(PROT_REL = PROT)
		if len(gene_feat_list) == 0:
			return
		
		feat_dict = {}
		feat_dict['gag'] = (colors.blue, 1)
		feat_dict['pol'] = (colors.red, -1)
		feat_dict['vif'] = (colors.green, 1)
		feat_dict['vpr'] = (colors.chocolate, -1)
		feat_dict['tat'] = (colors.orange, 1)
		feat_dict['rev'] = (colors.magenta, -1)
		feat_dict['env'] = (colors.darkorchid, 1)
		feat_dict['nef'] = (colors.pink, -1)
		feat_dict['vpu'] = (colors.purple, 1)
		
		
		for feat in gene_feat_list:
			feat.strand = feat_dict[feat.id][1]
			INPUT_TRACK.add_feature(feat, colour = feat_dict[feat.id][0])
			
	
	def RenderDiagram(self):
		"""
		RederDiagram
			A set of functions for annotating all of the features that have
			been annotated to the genome.  Overwrites any previous rendering.
		"""
		self.this_genome = GDDiagram(self.seq_name)
		
		wanted_tracks = ['Genes', 'HomIsland']
		track_dict = {}
		
		for i in range(len(wanted_tracks)):
			print str(i) + wanted_tracks[i]
			a_new = self.this_genome.new_track(i + 1, name = wanted_tracks[i],  
						greytrack = 1).new_set('feature')
			track_dict[wanted_tracks[i]] = a_new
		
		self.WriteGenesToDiagram(track_dict['Genes'])
		
		for this_annot in self.feature_annot:
			this_type = this_annot.type
			track_dict[this_type].add_feature(this_annot.GetSeqFeature())
			
		self.this_genome.draw()
	
	def LogAnnotation(self, TYPE, POS, SEQ, HOM, NAME):
		"""
		LogAnnotation
			Logs various types of annotations about the sequence for various
			features such as Genes, homology, and ELMs ... use none if the data
			does not apply
		"""
		
		if (type(POS) != types.TupleType):
			try:
				type(POS.type) == types.StringType
			except AttributeError:
				raise TypeError, 'If POS is not a tuple then it must be an ANNOT'
			new_feat = POS
		elif TYPE == 'Annot':
			new_feat = Annot(NAME, POS[0], POS[1], None)
		elif TYPE == 'HomIsland':
			new_feat = HomIsland(SEQ, POS[0], POS[1], HOM)
		elif TYPE == 'MIRNA':
			new_feat = HumanMiRNA(NAME, POS[0], POS[1], HOM)
		elif TYPE == 'ELM':
			new_feat = ELM(NAME, POS[0], POS[1], None)
			self.feature_annot_type['ELM'] = True
			logging.debug('ELM %(elm)s found on %(seq)s' % {'elm':str(new_feat), 'seq':self.seq_name})
		elif TYPE == 'TF':
			new_feat = TFSite(NAME, POS[0], POS[1], HOM)
		elif TYPE == 'ResistSite':
			new_feat = ResistSite(NAME, POS[0], POS[1], HOM)
		else:
			raise KeyError
		
		bisect.insort(self.feature_annot, new_feat)
		self.feature_annot_type[new_feat.type] = True
		
		
	def FindEqAnnot(self, WANTED_ANNOT, POS_FUDGE, IS_LIST = False):
		"""
		Attemps to find the equivelant Annotation recored within 
		self.feature_annot
		"""
		
		if IS_LIST:
			output = []
			for this_annot in WANTED_ANNOT:
				val = self.FindEqAnnot(this_annot, POS_FUDGE)
				output.append(val)
			return output
					
					
					
		type_fil = lambda x: x.type == WANTED_ANNOT.type
		name_fil = lambda x: x.name == WANTED_ANNOT.name
		pos_fun = lambda x: abs(x.start - WANTED_ANNOT.start)
		
		fil_feat_1 = filter(type_fil, self.feature_annot)
		fil_feat_2 = filter(name_fil, fil_feat_1)
		
		if len(fil_feat_2) > 0:
			pos_val = min(fil_feat_2, key = pos_fun)
		else:
			return
		
		if pos_fun(pos_val) < POS_FUDGE:
			return pos_val
		else:
			return
		
		
	def HumanMiRNAsite(self, CALIB_DICT, P_VAL = 0.05, NUM_THREADS = 10, ONLY_CHECK = None):
		"""
		HumanMiRNAsite
			Annotates the sites where human miRNA could potentially bind.  
			Makes calls to LogAnnotation with binding sites for latter
			descriptions.
		"""
		
		def HybridWorker(MI_RNA_NAME, THIS_CALIB, SEQ):
			"""
			HybridWorker
				A worker function which calls the HybridSeq function in a 
				threaded manner. Aquires lock a before making calls to 
				LogAnnotation to make sure there is no difficulty with
				the dictionary in this_ref.
			"""
			output_list = HybridSeq(THIS_CALIB[0], THIS_CALIB[1:], SEQ)
			if len(output_list) == 0:
				worker_seph.release()
				return
			with annot_lock:
				#need a lock because we are modifying lists and dictionaries 
				#in weird ways on the other side of the function call
				for this_spot in output_list:
					if this_spot[2] > P_VAL:
						continue
					this_val = (this_spot[0], 
									this_spot[0] + len(THIS_CALIB[0]))
					
					self.LogAnnotation('MIRNA', this_val, None, this_spot[2], 
										MI_RNA_NAME)
			worker_seph.release()
			return
			
		RNA_HYBRID_PATH = 'C:\\RNAHybrid\\'
		WANTED_THREADS = NUM_THREADS
		
		annot_lock = threading.Lock()
		worker_seph = threading.Semaphore(WANTED_THREADS)
		all_threads = []
		
		if ONLY_CHECK != None:
			m_rna_frac = ONLY_CHECK
			logging.debug('Only checking %(num)d miRNAs' % {'num':len(m_rna_frac)})
		else:
			m_rna_frac = CALIB_DICT.keys()
			logging.debug('Checking all %(num)d miRNAs' % {'num':len(m_rna_frac)})
		
		for this_mi_rna in m_rna_frac:
			worker_seph.acquire()
			this_thread = threading.Thread(target = HybridWorker, 
									args = (this_mi_rna, 
											CALIB_DICT[this_mi_rna], 
											self.my_sequence))
			this_thread.start()
			all_threads.append(this_thread)
		#make sure all threads have finished before exiting
		for this_thread in all_threads:
			this_thread.join()
		#make sure is set to true otherwise when checking anything that has 
		#NO miRNAs will be re-checked because it was never logged
		self.feature_annot_type['MIRNA'] = True
	def FindELMs(self, ELM_DICT):
		"""
		Finds the ELMs in the PROTIEN sequence and then using LogAnnotation 
		to mark thier locations on the NUCLEOTIDE sequence.
		"""
		for this_gene_name in self.annotation:
			this_gene = self.annotation[this_gene_name]
			logging.debug('%(seq)s has %(gene)s' % {'seq':self.seq_name, 'gene':this_gene.aa_seq})
			start_pos = this_gene.start
			rel_start = this_gene.rel_start
			for this_ELM in ELM_DICT:
				spot = ELM_DICT[this_ELM][1].search(this_gene.aa_seq)
				while spot != None:
					
					new_elm = ELM(this_ELM, start_pos + spot.start(), 
									start_pos + spot.end(), None)
					new_elm.GeneAnnot(this_gene_name, 
									spot.start(), spot.end())
					
					self.LogAnnotation('ELM', new_elm, None, None, this_ELM)
					
					spot = ELM_DICT[this_ELM][1].search(this_gene.aa_seq, 
														spot.start() + 1)
		self.feature_annot_type['ELM'] = True
	def FindHomIslands(self, REF_SEQ):
		"""
		Finds the homology islands that are identical to the reference 
		provided and then calls LogAnnotation.
		"""
		check_fun = lambda x: x.type == 'HomIsland'
		
		for this_annot in filter(check_fun, REF_SEQ.feature_annot):
			spot = string.find(self.my_sequence, this_annot.seq)
			while spot != -1:
				self.LogAnnotation('HomIsland', (spot, spot + len(this_annot.seq)),
									this_annot.seq, this_annot.hom, None)
				spot = string.find(self.my_sequence, this_annot.seq, spot+1)
		self.feature_annot_type['HomIsland'] = True
	
	def FindTFSites(self):
		"""
		Makes a call to AnnotUtils.TFChecker to annotate the transcription 
		factor binding sites.  It will then call LogAnnotation to add the data
		to the object.
		"""
		output = TFChecker(self.my_sequence)
		
		if len(output) > 0:
			for this_annot in output:
				spot = (this_annot[0], this_annot[0] + len(this_annot[4]))
				self.LogAnnotation('TF', spot, this_annot[4], 
									this_annot[3], this_annot[6])
		self.feature_annot_type['TF'] = True
	
	def FindPFAMSites(self, DICT_OBJ, DICT_LOCK):
		"""
		Makes a call to AnnotUtils.PFAMChecker to annotate the transcription 
		factor binding sites.  It will then call LogAnnotation to add the data
		to the object.
		"""
		output = PFAMChecker(self.annotation.values(), DICT_OBJ, DICT_LOCK)
		
		if len(output) > 0:
			for this_annot in output:
				self.LogAnnotation('PFAM', this_annot, None, None, None)
		self.feature_annot_type['PFAM'] = True
	
	def FindResistSites(self, RESIST_DICT):
		"""
		Finds the Resist Sites in the PROTIEN sequence and then using 
		LogAnnotation to mark thier locations on the NUCLEOTIDE sequence.
		"""
		
		for this_gene_name in self.annotation:
			this_gene = self.annotation[this_gene_name]
			logging.debug('%(seq)s has %(gene)s' % {'seq':self.seq_name, 'gene':this_gene.aa_seq})
			start_pos = this_gene.start
			rel_start = this_gene.rel_start
			for this_resist in RESIST_DICT:
				# if this_gene_name != this_resist[0]:
					# continue
				spot_wt = RESIST_DICT[this_resist][1].search(this_gene.aa_seq)
				spot_mut = RESIST_DICT[this_resist][2].search(this_gene.aa_seq)
				if spot_wt != None:
					new_resist = ResistSite(this_resist, start_pos + spot_wt.start(), 
									start_pos + spot_wt.end(), True)
					new_resist.GeneAnnot(this_gene_name, 
									spot_wt.start(), spot_wt.end())
					logging.debug('DRUG found WT spot')
				elif spot_mut != None:
					new_resist = ResistSite(this_resist, start_pos + spot_mut.start(), 
									start_pos + spot_mut.end(), False)
					new_resist.GeneAnnot(this_gene_name, 
									spot_mut.start(), spot_mut.end())
					logging.debug('DRUG found MUT spot')
				else:
					continue
				
				self.LogAnnotation('ResistSite', new_resist, None, None, 
										this_resist)
				
				
		self.feature_annot_type['ResistSite'] = True
			
	
	def CheckFeatures(self, FEAT_LIST, FUDGE_FACTOR, NUM_CHECKS = 3, 
									CHECK_CUTOFF = 1):
		"""
		Checks a list of provided features and returns a boolean array 
		indicating whether that feature is present on this sequence.
		"""
		
		all_checks = numpy.zeros((len(FEAT_LIST), NUM_CHECKS))
		
		min_func = lambda x: x[0]
		
		for this_check_num in xrange(NUM_CHECKS):
			#do multiple checks with different permutations
			#it may be slightly order dependent
			check_feat = copy.deepcopy(self.feature_annot)
			numpy.random.shuffle(check_feat)
			
			#loop through each feature in FEAT_LIST
			for this_feat in xrange(len(FEAT_LIST)):
				poss_items = []
				#loop through each feature in the shuffled copy of this list
				for this_check in xrange(len(check_feat)):
					val = check_feat[this_check].FuzzyEquals(FEAT_LIST[this_feat],
													FUDGE_POS = FUDGE_FACTOR)
					
					#val == -1 refers to "wrong type or not within "
					if val != -1:
						logging.debug('Checking %(selfname)s %(this)s:%(checkname)s %(test)s' \
							% {'this':str(check_feat[this_check]), 'test':str(FEAT_LIST[this_feat]),
								'selfname':self.seq_name, 'checkname':'TEST'})
						logging.debug('MARKED AS CORRECT')
						poss_items.append((val, this_check))
					else:
						pass
						#logging.debug('DISCARDED')
				
				if len(poss_items) > 0:
					best_ind = min(poss_items, key = min_func)[1]
					all_checks[this_feat, this_check_num] = 1
					check_feat.pop(best_ind)
					
				
		
		avg_find = numpy.sum(all_checks, axis=1)
		
		return avg_find > CHECK_CUTOFF
		
	def WriteFeatures(self, FEAT_LIST, FILE_HANDLE, FUDGE_FACTOR = 50, 
						NUM_CHECKS = 3, CHECK_CUTOFF = 1):
		"""
		Writes the feature list to a file-handle
		"""
		
		found_feats = self.CheckFeatures(FEAT_LIST, FUDGE_FACTOR, 
											NUM_CHECKS = NUM_CHECKS,
											CHECK_CUTOFF = CHECK_CUTOFF)
		
		trans_dict = {True: '1', False: '0'}
		
		found_feats = map(lambda x: trans_dict[x], found_feats.tolist())
		
		FILE_HANDLE.write(string.join(found_feats, '\t') + '\n')
		
	def DetermineResponder(self, METHOD):
		"""
		A stub function to provide compatibility with a PatSeq.
		ALWAYS RETURNS 2!!!!!
		"""
		return 2
		
	
	
class BkgSeq(ViralSeq):
	"""
	Used for Background Sequences
	"""

class RefSeq(ViralSeq):
	def __init__(self, FILENAME):
		self.my_sequence = None
		self.tested_subtype = None
		self.seq_name = None
		self.annotation = None
		#a windowed homology over the entire genome
		self.global_hom = None
		#identification of specific homology islands
		self.feature_annot = []
		self.multi_feature_annot = {}
		
		self.this_genome = None
		self.multi_genome = None
		
		self.feature_annot = []
		self.feature_annot_type = {}
		
		self.ParseGenbank(FILENAME)
		#self.DetSubtype()


	def ParseGenbank(self, FILENAME):
		with open(FILENAME, mode='r') as handle:
			this_record = SeqIO.parse(handle, 'genbank').next()

		self.my_sequence = this_record.seq.tostring()
		self.seq_name = this_record.id
		aa_start = 0
		self.annotation = {}
		for feat in this_record.features:
			if feat.qualifiers.has_key('gene') & feat.qualifiers.has_key('translation'):
				gene_name = feat.qualifiers['gene'][0]
				if '-' in gene_name:
					continue
				start_pos = feat.location._start.position
				end_pos = feat.location._end.position
				trans_data = feat.qualifiers['translation'][0]
				nt_seq = self.my_sequence[start_pos:end_pos]
				
				prod_name = feat.qualifiers['product'][0]
				
				this_gene = Gene(gene_name, prod_name, nt_seq,
								 trans_data, start_pos, end_pos)
				this_gene.GeneAnnot(None, aa_start, aa_start + len(trans_data))
				aa_start += len(trans_data)
				self.annotation[feat.qualifiers['gene'][0]] = this_gene
				
		
	def AnnotGlobalHom(self):
		"""
		AnnotGlobalHom
			Uses the data in self.global_hom to create a continious value of 
			homology
		"""
		if self.this_genome == None:
			self.this_genome = GDDiagram(self.seq_name)
		
		hom_track = self.this_genome.new_track(4, greytrack = 1, name = 'Homology')
		hom_set = hom_track.new_set('graph')
		hom_set.new_graph(self.global_hom, 'Homology', style = 'line')
		
		self.this_genome.draw()
		
	def DrawMultiGenome(self):
		
		self.multi_genome = GenomeDiagram(self.seq_name)
		
		for this_annot_list in self.multi_feature_annot.values():
			this_track = self.multi_genome.new_track(1)
			this_set = this_track.new_set('feature')
			for this_annot in this_annot_list:
				this_set.add_feature(this_annot.GetSeqFeature(), color = this_annot.color)
	
	def GetHomIslandFromFile(self, F_HANDLE):
		"""
		Read previously found Homology Islands from a file.
		"""
		
		island_dict = {}
		for this_line in F_HANDLE:
			island_dict[string.split(this_line, '/t')[0]] = None
		
		for this_island in island_dict:
			pos = self.my_sequence.find(this_island)
			if pos != -1:
				self.LogAnnotation('HomIsland', (pos, pos+len(this_island)), 
									this_island, 1.0, None)

class RefBase():
	def __init__(self, SOURCE_DIR, DEST_DIR, BUILD = False,
				 BLAST_DIR = 'C:\\local_blast\\', SEQ_FILES = None):


		self.ref_seqs = []
		if SOURCE_DIR != None:
			for this_file in filter(lambda x: x[-3:] == '.gb',
									os.listdir(SOURCE_DIR)):
					
				self.ref_seqs.append(RefSeq(SOURCE_DIR + this_file))
		else:
			for this_file in SEQ_FILES:
				self.ref_seqs.append(RefSeq(this_file))
		
		self.nt_name = DEST_DIR + 'ref_nt.fasta'
		self.aa_name = DEST_DIR + 'ref_aa.fasta'
		self.dir_name = DEST_DIR
		self.blast_dir = BLAST_DIR

		if BUILD:
			self.BuildDatabase()

	def __iter__(self):
		return iter(self.ref_seqs)
	
	def __getitem__(self, KEY):
		return self.GetRefSeq(KEY)
	
	def GetRefSeq(self, WANTED_REF):
		if WANTED_REF == None:
			return None
		for this_ref in self.ref_seqs:
			if this_ref.seq_name == WANTED_REF:
				return this_ref
		raise KeyError, "Couldn't find the RefSeq %(k)s" % {'k':WANTED_REF} 
	
	def FinalizeAnnotations(self):
		"""
		FinalizeAnnotations
			Occasionally the genbank parser is not able to read all of the gene
			sections of the genome.  FinalizeAnnotations uses the known
			annotations to infer the missing annotations and re-build the BLAST
			databases.
		"""
		
		for this_ref in self.ref_seqs:
			this_ref.TranslateAll(self)
			
		self.BuildDatabase()
	
	def BuildDatabase(self):
		"""
		BuildDatabase
				Uses the RefSeqs provided to create a protein and nucleotide
		blast database.
		"""
		logging.debug('Building Database')
		this_iter = IT.imap(RefSeq.GenomeToBioPython, iter(self.ref_seqs))
		
		logging.debug('Writting nt_file:' + self.nt_name)
		with open(self.nt_name, mode='w') as handle:
			SeqIO.write(this_iter, handle, "fasta")
		
		logging.debug('Writting aa_file:' + self.aa_name)
		with open(self.aa_name, mode='w') as handle:
			for this_data in self.ref_seqs:
				this_iter = this_data.AnnotToBioPython()
				SeqIO.write(this_iter, handle, "fasta")
		
		logging.debug('Making BLAST nt database:' + self.nt_name)
		command = self.MakeBLASTCommand('formatdb',
										self.nt_name,
										EXTRA_FLAGS = ' -p F -o T')
		
		ProcessVar = subprocess.Popen(command, shell = True)
		ProcessVar.wait()
		
		logging.debug('Making BLAST aa database:' + self.aa_name)
		command = self.MakeBLASTCommand('formatdb',
										self.aa_name,
										EXTRA_FLAGS = ' -p T -o T')

		ProcessVar = subprocess.Popen(command, shell = True)
		ProcessVar.wait()
		
		logging.debug('BLAST databases properly created')
		

	def BLASTx(self, INPUT_SEQ_RECORD, NUM_THREADS = None):
		"""
		BLASTx
				Performs a BLASTx query on the database to determine
		the possible translations of the INPUT_GENOME
		"""
		if NUM_THREADS == None:
			b_controler = PyBLAST.BLASTController(INPUT_SEQ_RECORD, 'BLASTx', 
										NUM_SEQS = 1, 
										BLAST_DIR = self.blast_dir,
										SCRATCH_DIR = self.dir_name)
			return b_controler.GetSingleResult()
		else:
			b_controler = PyBLAST.BLASTController(INPUT_SEQ_RECORD, 'BLASTx', 
										BLAST_DIR = self.blast_dir,
										SCRATCH_DIR = self.dir_name, 
										NUM_THREADS = NUM_THREADS)
			b_controler.start()
			return b_controler.ResultGen()

	def BLASTn(self, INPUT_SEQ_RECORD, NUM_THREADS = None):
		"""
		BLASTn
				Performs a BLASTn query to align the provided genome to
		the reference genomes.
		"""
		if NUM_THREADS == None:
			b_controler = PyBLAST.BLASTController(INPUT_SEQ_RECORD, 'BLASTn', 
										NUM_SEQS = 1, 
										BLAST_DIR = self.blast_dir,
										SCRATCH_DIR = self.dir_name)
			return b_controler.GetSingleResult()
		else:
			b_controler = PyBLAST.BLASTController(INPUT_SEQ_RECORD, 'BLASTn', 
										BLAST_DIR = self.blast_dir,
										SCRATCH_DIR = self.dir_name, 
										NUM_THREADS = NUM_THREADS)
			b_controler.start()
			return b_controler.ResultGen()
			
	def MakeBLASTCommand(self, BLAST_TYPE, INPUT_FILE, EXTRA_FLAGS = None):
		"""
		MakeBLASTCommand
				Creates a command that can be passed to Popen for
		executing blast commands ... supports: formatdb, blastn,
		blastx.
		"""
		if BLAST_TYPE == 'formatdb':
			command = self.blast_dir + 'bin\\formatdb.exe'

		command += ' -i ' + INPUT_FILE
		
		if EXTRA_FLAGS != None:
			command +=  EXTRA_FLAGS

		return command



class BLASTContext:
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


