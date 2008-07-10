#-----------------------------------------
#  Author:     Perry Evans
#              evansjp@mail.med.upenn.edu
#  2008
#
#-----------------------------------------
"""
Functions for dealing with FASTA files.  The
format I'm using is:

>gene_name
seq
seq
seq
... .
"""

def printFASTA_forGenes(gene_file, fasta_file, output_file):
    """ Given a file of genes and a FASTA file,
        print the fasta entries in the gene file.
    
    @param gene_file: one gene per line; you want seqs for these genes
    @param fasta_file: seqs
    @param output_file: where to stick the seqs
    """
    
    genes = {}
    gene_f = open(gene_file)
    for line in gene_f:
        genes[line.strip()] = True
    gene_f.close()

    fasta_f = open(fasta_file)
    fout = open(output_file, 'w')
    line = fasta_f.readline()
    while line != '':
        if line[0] == '>':
            if genes.has_key(line[1:].strip()):
                fout.write(line)
                line = fasta_f.readline()
                while line[0] != '>':
                    fout.write(line)
                    line = fasta_f.readline()
                    if line == '':
                        break
            else:
                line = fasta_f.readline()
                while line[0] != '>':
                    line = fasta_f.readline()
                    if line == '':
                        break
    fasta_f.close()
    fout.close()

def loadFASTA(fasta_file):
    """ Make a {} from a FASTA file.

    @param fasta_file: >name
                       seq
                       seq ...
    @return: fasta[gene] = seq (oneline)
    """

    fasta = {}
    fasta_f = open(fasta_file)
    name = ''
    seq = ''
    line = fasta_f.readline()
    while line != '':
        if line[0] == '>':
            if name != '':
                fasta[ name ] = seq
            name = line[1:].strip()
            seq = ''
        else:
            seq = seq + line.strip()
        line = fasta_f.readline()
    fasta_f.close()

    if name != '':    
        fasta[name] = seq

    return fasta

def prettyPrint_runner(geneid, seq_string):
    """ Return FASTA formatted version of 
        this gene and its sequence.

    @param geneid: gene_name
    @param seq_string: FASTA seq on one line w/ no whitespace
    @return: formatted FASTA version
    """

    breakcount = 50
    line = '>' + geneid + '\n'
    while len(seq_string) > breakcount:
        line = line + seq_string[0:breakcount] + '\n'
        seq_string = seq_string[breakcount:]
    if seq_string != '':
        line = line + seq_string
    return line

def prettyPrint(geneid, seq_string):
    """ Print FASTA formatted version of this
        gene and its sequence.

    @param geneid: gene_name
    @param seq_string: FASTA seq on one line w/ no whitespace
    """

    print prettyPrint_runner(geneid, seq_string)

def prettyPrint_file(geneid, seq_string, afile):
    """ Write FASTA formatted gene/seq to this open file.

    @param geneid: gene_name
    @param seq_string: FASTA seq on one line w/ no whitespace
    @param afile: an open output file
    """
    
    afile.write( prettyPrint_runner(geneid, seq_string) + '\n')
