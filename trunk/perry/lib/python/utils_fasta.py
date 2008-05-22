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

def loadFASTA(fasta_file):
    fasta = {}
    f = open(fasta_file)
    name = ''
    line = f.readline()
    while line != '':
        if line[0] == '>':
            if name != '':
                fasta[ name ] = seq
            name = line[1:].strip()
            seq = ''
        else:
            seq = seq + line.strip()
        line = f.readline()
    f.close()

    if name != '':    
        fasta[name ] = seq

    return fasta

def prettyPrint_runner(geneid, seq_string):
    breakcount = 50
    line = '>' + geneid + '\n'
    while len(seq_string) > breakcount:
        line = line + seq_string[0:breakcount] + '\n'
        seq_string = seq_string[breakcount:]
    if seq_string != '':
        line = line + seq_string
    return line

def prettyPrint(geneid, seq_string):
    print prettyPrint_runner(geneid, seq_string)

def prettyPrint_file(geneid, seq_string, file):
    file.write( prettyPrint_runner(geneid, seq_string) + '\n')
