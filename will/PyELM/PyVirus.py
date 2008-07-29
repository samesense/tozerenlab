from __future__ import with_statement
import itertools as IT
import HIVGenoTypeChecker as GenoTyping
import subprocess
import tempfile
import threading
import time
import os
import re
import types
import Queue
import PyBLAST
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

class Annot():
	def __init__(self, NAME, START_POS, END_POS):
		self.name = NAME
		self.start = START_POS
		self.end = END_POS
		self.type = 'default'
	
	def GetSeqFeature(self):
		"""
		Returns a BioPython SeqFeature describing the gene represented
		"""
		seq_loc = SeqFeatClass.FeatureLocation(self.start, self.end)
		seq_feat = SeqFeature(location = seq_loc, strand = 1, 
								id = self.name)
		return seq_feat

class Gene(Annot):
	def __init__(self, NAME, PRODUCT, NT_SEQ, AA_SEQ, START, END):
		self.nt_seq = NT_SEQ
		self.aa_seq = AA_SEQ
		self.start = START
		self.end = END
		self.name = NAME
		self.product = PRODUCT
		self.type = 'GENE'

	def __str__(self):
		temp_str = '(Gene Class '
		temp_str += 'Name: ' + self.name
		temp_str += ', Start: ' + str(self.start) + ')'
		return temp_str
		
	def __hash__(self):
		return hash(self.aa_seq)
		
class HumanMiRNA(Annot):
	def __init__(self, NAME, START_POS, END_POS):
		self.name = NAME
		self.start = START_POS
		self.end = END_POS
		self.type = 'HumanMiRNA'
		self.color = colors.blue
		
class HomIsland(Annot):
	def __init__(self, SEQ, START_POS, HOM):
		self.seq = SEQ
		self.start = START_POS
		self.hom = HOM
		self.name = SEQ
		self.end = START_POS + len(SEQ)
		self.type = 'HomIsland'
		self.color = colors.green
		
	def __hash__(self):
		return hash(self.seq)


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
		self.this_genome = None

	def __hash__(self):
		return hash(self.my_sequence + self.seq_name)

	def __getattr__(self, NAME):
		if NAME == 'tested_subtype':
			self.DetSubtype()
			return self.tested_subtype
		else:
			raise AttributeError
		
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
	
	def GetSeqFeatures(self):
		"""
		GetSeqFeatures
			Returns a list of SeqFeatures describing each protein in 
			self.annotation suitalbe for using with GenomeDiagram
		"""
		feat_list = []
		for this_gene in self.annotation:
			feat_list.append(self.annotation[this_gene].GetSeqFeature())
		
		return feat_list

	def DetSubtype(self):
		"""
		DetSubtype
				Determines the subtype of the sequence using the
		HIVSubtype module.
		"""
		self.tested_subtype = GenoTyping.GetSimple(self.my_sequence)

	def TranslateAll(self, REFBASE):
		"""
		TranslateAll
				Uses a BLASTx query to determine the which ORFs
		correspond to HIV-1 proteins.  It then stores the data in a
		dictionary keyed by "genename".
		"""
		blast_record = REFBASE.BLASTx(self.GenomeToBioPython())
		gene_dict = {}
		reg = re.compile('.*?:(\w{3}).*')
		for this_check in blast_record.alignments:
			this_gene = str(reg.match(this_check.title).groups()[0])
			if not(gene_dict.has_key(this_gene)):
				aa_seq = str(this_check.hsps[0].query)
				start_pos = this_check.hsps[0].query_start
				end_pos = start_pos + this_check.length*3
				nt_seq = self.my_sequence[start_pos:end_pos]
				gene = Gene(this_gene, this_gene.upper(),
							nt_seq, aa_seq, start_pos, end_pos)
				gene_dict[this_gene] = gene
		self.annotation = gene_dict
	
	def CheckSeq(self, INPUT_SEQ):
		"""
		CheckSeq
			Returns True if the sequence is present within the sequence and False
			otherwise.
		"""
		return self.my_sequence.find(INPUT_SEQ) != -1
	
	def WriteGenesToDiagram(self, INPUT_TRACK):
		"""
		WriteGenomeDiagram
			Use the GenomeDiagram module to write a publication quality figure
			of the genes and features in self.annotation
		"""
		
		gene_feat_list = self.GetSeqFeatures()
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
	
	def LogAnnotation(self, TYPE, DATA):
		"""
		LogAnnotation
			Logs various types of annotations about the sequence for various
			features such as Genes, homology, and ELMs
		"""
		
		if TYPE == 'Annot':
			self.feature_annot.append(Annot(DATA[0], DATA[1], DATA[2], DATA[4]))
		elif TYPE == 'HomIsland':
			self.feature_annot.append(HomIsland(DATA[0], DATA[1], DATA[2]))
		else:
			raise KeyError
		
		
		
class BkgSeq(ViralSeq):
	"""
	Used for Background Sequences
	"""


class PatSeq(ViralSeq):
	def __init__(self, TEST_SEQ, SEQ_NAME):
		self.my_sequence = TEST_SEQ
		self.pat_data = None
		self.tested_subtype = None
		self.seq_name = SEQ_NAME
		self.this_genome = None

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
		
		self.ParseGenbank(FILENAME)
		#self.DetSubtype()


	def ParseGenbank(self, FILENAME):
		with open(FILENAME, mode='r') as handle:
			this_record = SeqIO.parse(handle, 'genbank').next()

		self.my_sequence = this_record.seq.tostring()
		self.seq_name = this_record.id
		
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
	
	def AnnotMultiGenome(self, TEST_NAME, FEATURE_NAME, 
							FEATURE_POS, FEATURE_TYPE):
		"""
		AnnotMultiGenome
			Annotates a figure where each TEST mapping has a different cirlce
			and each feature type has a different color.  Helpful for seeing 
			patterns across multiple mappings in the positioning of features
		"""
		
		annot_list = self.multi_feature_annot.setdefault(TEST_NAME, [])
		
		if FEATURE_TYPE == 'MIRNA':
			this_annot = HumanMiRNA(FEATURE_NAME, FEATURE_POS[0], FEATURE_POS[1])
			annot_list.append(this_annot)
			self.multi_feature_annot[TEST_NAME] = annot_list
		else:
			raise KeyError
		
		
	def DrawMultiGenome(self):
		
		self.multi_genome = GenomeDiagram(self.seq_name)
		
		for this_annot_list in self.multi_feature_annot.values():
			this_track = self.multi_genome.new_track(1)
			this_set = this_track.new_set('feature')
			for this_annot in this_annot_list:
				this_set.add_feature(this_annot.GetSeqFeature(), color = this_annot.color)
		

class RefBase():
	def __init__(self, SOURCE_DIR, DEST_DIR, BUILD = False,
				 BLAST_DIR = 'C:\\local_blast\\'):


		self.ref_seqs = []
		for this_file in filter(lambda x: x[-3:] == '.gb',
								os.listdir(SOURCE_DIR)):
				
			self.ref_seqs.append(RefSeq(SOURCE_DIR + this_file))
		
		self.nt_name = DEST_DIR + 'ref_nt.fasta'
		self.aa_name = DEST_DIR + 'ref_aa.fasta'
		self.dir_name = DEST_DIR
		self.blast_dir = BLAST_DIR

		if BUILD:
			self.BuildDatabase()

	def __iter__(self):
		return iter(self.ref_seqs)
	
	def GetRefSeq(self, WANTED_REF):
		for this_ref in self.ref_seqs:
			if this_ref.seq_name == WANTED_REF:
				return this_ref
		raise KeyError
	
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

		this_iter = IT.imap(RefSeq.GenomeToBioPython, iter(self.ref_seqs))
		with open(self.nt_name, mode='w') as handle:
			SeqIO.write(this_iter, handle, "fasta")

		with open(self.aa_name, mode='w') as handle:
			for this_data in self.ref_seqs:
				this_iter = this_data.AnnotToBioPython()
				SeqIO.write(this_iter, handle, "fasta")

		command = self.MakeBLASTCommand('formatdb',
										self.nt_name,
										EXTRA_FLAGS = ' -p F -o T')
		
		ProcessVar = subprocess.Popen(command, shell = True)
		ProcessVar.wait()

		command = self.MakeBLASTCommand('formatdb',
										self.aa_name,
										EXTRA_FLAGS = ' -p T -o T')

		ProcessVar = subprocess.Popen(command, shell = True)
		ProcessVar.wait()

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


