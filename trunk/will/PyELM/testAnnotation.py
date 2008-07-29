from __future__ import with_statement
import nose.tools
import HIVDatabase
import PyVirus
import os
import time
import numpy
from Bio import SeqIO
import itertools
import re

def setupModule():
	"""
	Determine if the required files are present ... otherwise create 
	them.
	"""
	
	dest_dir = os.environ['PYTHONSCRATCH']
	source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'
	
	PyVirus.RefBase(source_dir, dest_dir, BUILD = True)
	
	needed_files = ['test_self_KEEP_map.slf', 
						'test_self_KEEP_map.slf']
	
	present_files = os.listdir(dest_dir)
	
	for this_file in needed_files:
		if not(this_file in present_files):
			mapping_base = HIVDatabase.MappingBase(source_dir, dest_dir,
												'test_self_KEEP')
			with open(seq_file) as handle:
				seq_iter = SeqIO.parse(handle, 'fasta')
				mapping_base.AddtoShelf(seq_iter)
			break
	
	mirna_dict_file = 'RNAiCalibrations_KEEP.pkl'
	
	if not(mirna_dict_file in present_files):
		DIREC = os.environ['MYDOCPATH'] + 'HIVRNAi\\'
		DATABASENAME = 'human_predictions_4dwnload.txt'
		HIVDatabase.CalibrateRNAi(DIREC + DATABASENAME, 
									dest_dir + mirna_dict_file)

def testAnnotation():
	"""
	Test Annotations
	"""
	dest_dir = os.environ['PYTHONSCRATCH']
	source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

	mapping_base = HIVDatabase.MappingBase(source_dir, dest_dir,
									   'test_self_KEEP')
	mapping_base.ref_base.FinalizeAnnotations()
	this_iter = itertools.izip(iter(mapping_base.ref_base),
					itertools.repeat(None,5))
	for this_ref in this_iter:
		#yield CheckAnnotGlobalHom, mapping_base, this_ref[0]
		#yield CheckFindWindows, mapping_base, this_ref[0]
		#yield CheckMakeDiagram, this_ref[0]
		yield CheckHumanMiRNA, mapping_base, this_ref[0]
		yield CheckDrawMulti, this_ref[0]
	
def CheckAnnotGlobalHom(MAPPING_BASE, THIS_REF):
	
	dest_dir = os.environ['PYTHONSCRATCH'] + THIS_REF.seq_name + '_KEEP.pdf'
	MAPPING_BASE.WindowedHomology(THIS_REF.seq_name)
	
	nose.tools.assert_true(len(THIS_REF.global_hom) != 0,
				'Did not make Homology properly in:' + THIS_REF.seq_name)

def CheckFindWindows(MAPPING_BASE, THIS_REF):
	"""
	Test the FindWindows
	"""
	MAPPING_BASE.FindWindows(THIS_REF.seq_name, 15, .6)
	num_win = len(THIS_REF.feature_annot)
	
	nose.tools.assert_true(num_win > 0, 
			'No homology windows were found in: ' + THIS_REF.seq_name)

def CheckMakeDiagram(THIS_REF):
	
	THIS_REF.RenderDiagram()
	THIS_REF.AnnotGlobalHom()
	
	dest_dir = os.environ['PYTHONSCRATCH']
	dest_file = dest_dir + THIS_REF.seq_name + '_KEEP.pdf'
	
	simple_name = THIS_REF.seq_name + '_KEEP.pdf'
	
	THIS_REF.this_genome.write(dest_file, 'PDF')
	
	nose.tools.assert_true(simple_name in os.listdir(dest_dir), 
			'Did not make ' + THIS_REF.seq_name + '_KEEP.pdf') 

def CheckHumanMiRNA(MAPPING_BASE, THIS_REF):
	
	dest_dir = os.environ['PYTHONSCRATCH']
	mirna_dict_file = 'RNAiCalibrations_KEEP.pkl'
	
	MAPPING_BASE.HumanMiRNAsite(dest_dir + mirna_dict_file, THIS_REF.seq_name)
	
	num_feat = len(THIS_REF.multi_feature_annot)
	
	nose.tools.assert_true(num_feat > 0, 'Did not create any features.')

def CheckDrawMulti(THIS_REF):
	
	dest_dir = os.environ['PYTHONSCRATCH']
	file_name = THIS_REF.seq_name + '_multi_KEEP.pdf'
	THIS_REF.DrawMultiGenome()
	THIS_REF.multi_genome.write(dest_dir + file_name, 'PDF')
	
	
def tearDownModule():
	checker = re.compile('.*KEEP.*')
	time.sleep(1)
	file_list = os.listdir(os.environ['PYTHONSCRATCH'])
	for this_file in file_list:
		if checker.match(this_file) == None:
			os.remove(os.environ['PYTHONSCRATCH'] + this_file)

    