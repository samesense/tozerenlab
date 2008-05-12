# given a list of genes
# and a FASTA file,
# print the fasta entries
def printFASTA_forGenes(geneLsFile, FASTA_file, output_file):
    genes = {}
    f = open(genesFile)
    for line in f.xreadlines():
        genes[line.strip()] = True
    f.close()

    f = open(FASTA_file)
    fout = open(output_file, 'w')
    line = f.readline()
    while line != '':
        if line[0] == '>':
            if have.has_key(line[1:].strip()):
                fout.write(line)
                line = f.readline()
                while line[0] != '>':
                    fout.write(line)
                    line = f.readline()
                    if line == '': break
            else:
                line = f.readline()
                while line[0] != '>':
                    line = f.readline()
                    if line == '': break
    f.close()
    fout.close()
