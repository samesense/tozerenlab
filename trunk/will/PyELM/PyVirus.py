from __future__ import with_statement
import shelve
#import PyAlign
import itertools as IT
import HIVGenoTypeChecker as GenoTyping
import subprocess
import tempfile
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein, generic_nucleotide
from Bio.Blast import NCBIXML





def MakeMapping(ALIGNMENT):
	REF_ITER = iter(ALIGNMENT[1])
	TEST_ITER = iter(ALIGNMENT[0])
	filterFun = lambda x: x[0] != '-'
	mapping_val = []
	num_val = []
	ref_num_inter = IT.chain(IT.ifilter(filterFun, IT.izip(REF_ITER, IT.count())), IT.repeat(('-',len(ALIGNMENT[1]))))
	test_iter = IT.ifilter(filterFun, IT.izip(TEST_ITER, ref_num_inter))
	for this_iter in test_iter:
		mapping_val.append(this_iter[0])
		num_val.append(this_iter[1][1])
	print ALIGNMENT[1] + str(len(ALIGNMENT[1]))
	print ALIGNMENT[0]
	print str(mapping_val) + str(len(mapping_val))
	print str(num_val) + str(len(num_val))

class Gene():
        def __init__(self, NAME, PRODUCT, NT_SEQ, AA_SEQ, START, END):
                self.nt_seq = NT_SEQ
                self.aa_seq = AA_SEQ
                self.start = START
                self.end = END
                self.name = NAME
                self.product = PRODUCT

        def __str__(self):
                temp_str = '(Gene Class '
                temp_str += 'Name: ' + self.name
                temp_str += ', Start: ' + str(self.start) + ')'
                return temp_str


class ViralSeq():
        def __init__(self, TEST_SEQ, SEQ_NAME):
                self.my_sequence = TEST_SEQ
                self.tested_subtype = None
                self.seq_name = SEQ_NAME
                self.annotation = None

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

class BkgSeq(ViralSeq):
        def __init__(self, TEST_SEQ, SEQ_NAME):
                self.my_sequence = TEST_SEQ
                self.tested_subtype = None
                self.seq_name = SEQ_NAME


class PatSeq(ViralSeq):
        def __init__(self, TEST_SEQ, SEQ_NAME):
                self.my_sequence = TEST_SEQ
                self.pat_data = None
                self.tested_subtype = None
                self.seq_name = SEQ_NAME


class RefSeq(ViralSeq):
        def __init__(self, FILENAME):
                self.my_sequence = None
                self.known_subtype = None
                self.seq_name = None
                self.annotation = None

                self.ParseGenbank(FILENAME)

                

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
	


class RefBase():
        def __init__(self, SOURCE_DIR, DEST_DIR, BUILD = True,
                     BLAST_DIR = 'C:\\local_blast\\'):


                self.ref_seqs = []
                for this_file in filter(lambda x: x[-3:] == '.gb',
                                        os.listdir(SOURCE_DIR)):
                        
                        self.ref_seqs.append(RefSeq(SOURCE_DIR + this_file))
                
                self.nt_name = DEST_DIR + 'ref_nt.fasta'
                self.aa_name = DEST_DIR + 'ref_aa.fasta'
                self.dir_name = DEST_DIR
                self.blast_dir = BLAST_DIR

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


                


                command = self.blast_dir + 'bin\\formatdb.exe'
                command += ' -i ' +  self.nt_name
                command += ' -p F -o T'

                ProcessVar = subprocess.Popen(command, shell = True)
                ProcessVar.wait()

                command = self.blast_dir + 'bin\\formatdb.exe'
                command += ' -i ' + self.aa_name
                command += ' -p T -o T'

                ProcessVar = subprocess.Popen(command, shell = True)
                ProcessVar.wait()
                


        def BLASTx(self, INPUT_GENOME):
                """
                BLASTx
                        Performs a BLASTx query on the database to determine
                the possible translations of the INPUT_GENOME
                """
                with BLASTContext(self.dir_name, 2) as file_names:
                        print file_names
                                

                        with open(file_names[0], mode = 'w') as handle:
                                SeqIO.write([INPUT_GENOME],
                                            handle, 'fasta')


                        command = self.blast_dir + 'bin\\blastall.exe'
                        command += ' -p blastx -m 7'
                        command += ' -d ' + self.aa_name
                        command += ' -o ' + file_names[1]
                        command += ' -i ' + file_names[0]

                        ProcessVar = subprocess.Popen(command)
                        ProcessVar.wait()
                       

                        with open(file_names[1], mode = 'r') as handle:
                                this_blast = NCBIXML.parse(handle).next()

                        

                return this_blast


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

                
                




        
        


    
