from __future__ import with_statement
import nose.tools
import nose
import HIVDatabase
import PyVirus
import os
import time
import numpy
from Bio import SeqIO
import itertools
import re

import logging



def Slowsetup():
	"""
	Determine if the required files are present ... otherwise create 
	them.
	"""
	
	dest_dir = os.environ['PYTHONSCRATCH']
	source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	seq_file = os.environ['MYDOCPATH'] + 'PyELM\\500_seqs.fasta'
	
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

									
									
def setupModule():
	"""
	Determine if the required files are present ... otherwise create 
	them.
	"""
	dest_dir = os.environ['PYTHONSCRATCH']
	source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	seq_file = os.environ['MYDOCPATH'] + 'PyELM\\500_seqs.fasta'
	
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

def tearDownModule():
	checker = re.compile('.*KEEP.*')
	time.sleep(1)
	file_list = os.listdir(os.environ['PYTHONSCRATCH'])
	for this_file in file_list:
		if checker.match(this_file) == None:
			os.remove(os.environ['PYTHONSCRATCH'] + this_file)

#@nose.with_setup(Fastersetup, tearDown)
def testAnnotation_FAST():
	"""
	Test Annotations
	"""
	dest_dir = os.environ['PYTHONSCRATCH']
	source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	seq_file = os.environ['MYDOCPATH'] + 'PyELM\\50_seqs.fasta'

	mapping_base = HIVDatabase.MappingBase(source_dir, dest_dir,
									   'test_self_KEEP')
	#mapping_base.ref_base.FinalizeAnnotations()
	this_iter = itertools.izip(iter(mapping_base.ref_base),
					itertools.repeat(None,1))
	for this_ref in this_iter:
		#yield CheckAnnotGlobalHom, mapping_base, this_ref[0]
		#yield CheckFindWindows, mapping_base, this_ref[0]
		#yield CheckMakeDiagram, this_ref[0]
		yield CheckDrawMulti, mapping_base, this_ref[0]

#@nose.with_setup(Slowsetup, tearDown)
def testAnnotation_SLOW():
	"""
	Test Annotations
	"""
	dest_dir = os.environ['PYTHONSCRATCH']
	source_dir = os.environ['MYDOCPATH'] + 'hivsnppredsvn\\HIVRefs\\'
	seq_file = os.environ['MYDOCPATH'] + 'PyELM\\500_seqs.fasta'

	mapping_base = HIVDatabase.MappingBase(source_dir, dest_dir,
									   'test_self_KEEP')
	#mapping_base.ref_base.FinalizeAnnotations()
	this_iter = itertools.izip(iter(mapping_base.ref_base),
					itertools.repeat(None,1))
	for this_ref in this_iter:
		#yield CheckAnnotGlobalHom, mapping_base, this_ref[0]
		#yield CheckFindWindows, mapping_base, this_ref[0]
		#yield CheckMakeDiagram, this_ref[0]
		yield CheckDrawMulti, mapping_base, this_ref[0]
	
def CheckAnnotGlobalHom(MAPPING_BASE, THIS_REF):
	
	dest_dir = os.environ['PYTHONSCRATCH'] + THIS_REF.seq_name + '_KEEP.pdf'
	MAPPING_BASE.WindowedHomology(THIS_REF.seq_name)
	
	nose.tools.assert_true(len(THIS_REF.global_hom) != 0,
				'Did not make Homology properly in:' + THIS_REF.seq_name)

def CheckFindWindows(MAPPING_BASE, THIS_REF):
	"""
	Test the FindWindows
	"""
	MAPPING_BASE.FindWindows(THIS_REF.seq_name, 15, .4)
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


def CheckDrawMulti(MAPPING_BASE, THIS_REF):
	
	dest_dir = os.environ['PYTHONSCRATCH']
	lin_file_name = THIS_REF.seq_name + '_lin_multi_KEEP.pdf'
	cir_file_name = THIS_REF.seq_name + '_cir_multi_KEEP.pdf'
	lin_prot_file_name = THIS_REF.seq_name + '_linPROT_multi_KEEP.pdf'
	cir_prot_file_name = THIS_REF.seq_name + '_cirPROT_multi_KEEP.pdf'
	
	if THIS_REF.seq_name == 'AF033819.3':
		#this_fig = MAPPING_BASE.MakeMultiDiagram(THIS_REF.seq_name, FORCE = True)
		this_fig = MAPPING_BASE.MakeMultiDiagram(THIS_REF.seq_name)
	else:
		this_fig = MAPPING_BASE.MakeMultiDiagram(THIS_REF.seq_name)
	
	
	this_fig[0].draw()
	this_fig[0].write(dest_dir + cir_file_name, 'PDF')
	
	this_fig[0].draw(format='linear', fragments = 1)
	this_fig[0].write(dest_dir + lin_file_name, 'PDF')
	
	this_fig[1].draw()
	this_fig[1].write(dest_dir + cir_prot_file_name, 'PDF')
	
	this_fig[1].draw(format='linear', fragments = 1)
	this_fig[1].write(dest_dir + lin_prot_file_name, 'PDF')
	

    