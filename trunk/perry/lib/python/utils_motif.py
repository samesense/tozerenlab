#---------------------------------------
#
# Author     Perry Evans
#            evansjp@mail.med.upenn.edu
# 2008
#---------------------------------------
"""
Functions for reading the my protein 
annotation file, which is organized:

geneid st stp annotation_name seq annotation_tool

and making annotated multiple alignment plots.
"""
import utils_graph
import pylab, math

def annotation2protein(afile, tool_dict):
    d = dict()
    f = open(afile)
    for line in f.xreadlines():
        [geneid, start, stop,
         domain, seq_desc, tool] = [x.strip() for x in line.split('\t')]
        start = int(start)
        stop = int(stop)
        if tool_dict.has_key( tool ):
            if not d.has_key( domain ): d[domain] = {}
            if not d[domain].has_key(geneid): d[domain][geneid] = []
            add = True
            for [st, stp, old_desc] in d[domain][geneid]:
                if start == st and stp == stop:
                    add = False
            if add:
                d[domain][geneid].append([start, stop, seq_desc])
    f.close()
    return d

def protein2annotation(afile, tool_dict):
    """ Parse the given file and return {}
        mapping proteins to their annotations.

    @param afile: annotation file; format
    geneid st stp annotation seq tool
    """

    d = dict()
    f = open(afile)
    for line in f.xreadlines():
        [geneid, start, stop,
         domain, seq_desc, tool] = [x.strip() for x in line.split('\t')]
        start = int(start)
        stop = int(stop)
        if tool_dict.has_key( tool ):
            if not d.has_key( geneid ): d[geneid] = {}
            if not d[geneid].has_key(domain):
                d[ geneid ][ domain ] = []
            add = True
            for [st, stp, old_desc] in d[geneid][domain]:
                if start == st and stp == stop:
                    add = False
            if add:
                d[geneid][domain].append([start, stop, seq_desc])
    f.close()
    return d

def protein2annotation_forMotifs(afile, motifs):
    d = dict()
    f = open(afile)
    for line in f.xreadlines():
        [geneid, start, stop,
         domain, seq_desc, tool] = [x.strip() for x in line.split('\t')]
        start = int(start)
        stop = int(stop)
        if motifs.has_key(domain):
            if not d.has_key( geneid ): d[ geneid ] = {}
            if not d[geneid].has_key(domain):
                d[ geneid ][ domain ] = []
            add = True
            for [st, stp, old_desc] in d[geneid][domain]:
                if start == st and stp == stop:
                    add = False
            if add:
                d[geneid][domain].append([start, stop, seq_desc])
    f.close()
    return d

def annotation2protein_forMotifs(afile, motifs):
    d = dict()
    f = open(afile)
    for line in f.xreadlines():
        [geneid, start, stop,
         domain, seq_desc, tool] = [x.strip() for x in line.split('\t')]
        start = int(start)
        stop = int(stop)
        if motifs.has_key(domain):
            if not d.has_key( domain ): d[ domain ] = {}
            if not d[domain].has_key(geneid):
                d[ domain ][ geneid ] = []
            add = True
            for [st, stp, old_desc] in d[domain][geneid]:
                if start == st and stp == stop:
                    add = False
            if add:
                d[domain][geneid].append([start, stop, seq_desc])
    f.close()
    return d

def motif2protein_forProteinLs(afile, ls):
    d = dict()
    f = open(afile)
    for line in f.xreadlines():
        sp = line.split('\t')
        if ls.has_key( sp[0] ):
            if not d.has_key( sp[3] ): d[ sp[3] ] = dict()
            d[ sp[3] ][ sp[0] ] = True
    f.close()
    return d

def mergeMotifs(start, end, ls):
    merged = False
    for i in xrange(len(ls)):
        oldStart = ls[i][0]
        oldEnd = ls[i][1]
        if start >= oldStart and end <= oldEnd:
            merged = True
        elif start <= oldStart and end >= oldEnd:
            ls[i] = [start, end]
            merged = True
        elif start <= oldStart and end >= oldStart and end <= oldEnd:
            ls[i] = [start, oldEnd]
            merged = True
        elif start >= oldStart and start <= oldEnd and end >= oldEnd:
            ls[i] = [oldStart, end]
            merged = True
    if not merged:
        ls.append( [start, end] )

def printAnnotation_runner(annotationLs, fasta):
    toColor = []
    for [start, stop, seq] in annotationLs:
        mergeMotifs(start, stop, toColor)
    line = ''
    toColor = sorted(toColor)    
    previousStop = 0
    for i in xrange( len(toColor) -1 ):        
        [start, stop] = toColor[i]
        [nextStart, nextStop] = toColor[i+1]
        line = line + fasta[previousStop:start-1] + '\033[1;31m' + fasta[start-1:stop] + '\033[0m'
        previousStop = stop
    [start, stop] = toColor[-1]
    line = line + fasta[previousStop:start-1] + '\033[1;31m' + fasta[start-1:stop] + '\033[0m' + fasta[stop:]
    return line

def printAnnotation(annotationLs, fasta):
    print printAnnotation_runner(annotationLs, fasta)

def printAnnotation_tofile(annotationLs, fasta, afile):
    afile.write(printAnnotation_runner(annotationLs, fasta) + '\n')

def printAnnotation01_runner(annotationLs, fasta):
    toColor = []
    for [start, stop, seq] in annotationLs:
        mergeMotifs(start, stop, toColor)
    line = ''
    toColor = sorted(toColor)    
    previousStop = 0
    for i in xrange( len(toColor) -1 ):        
        [start, stop] = toColor[i]
        [nextStart, nextStop] = toColor[i+1]
        for x in xrange(previousStop,start-1):
            line = line + '0 '
        for x in xrange(start-1,stop):
            line = line + '1 '
        #line = line + fasta[previousStop:start-1] + '\033[1;31m' + fasta[start-1:stop] + '\033[0m'
        previousStop = stop
    [start, stop] = toColor[-1]
    #line = line + fasta[previousStop:start-1] + '\033[1;31m' + fasta[start-1:stop] + '\033[0m' + fasta[stop:]
    for x in xrange(previousStop, start-1):
        line = line + '0 '
    for x in xrange(start-1, stop):
        line = line + '1 '
    for x in xrange(stop, len(fasta)):
        line = line + '0 '
    return line.strip()

def printAnnotation01(annotationLs, fasta):
    print printAnnotation01_runner(annotationLs, fasta)

def printAnnotation01_tofile(annotationLs, fasta, afile):
    afile.write(printAnnotation01_runner(annotationLs, fasta) + '\n')
        
def extendSeq(max_len, seq):
    """ All sequnces in the alignment must
        be the same length.  This fills in
        the seq with 0 until max_len is reached.
    
    @param max_len: desired seq len
    @param seq as a list of 1 and 0 for motif presence/absense
    @return: seq w/ filled in 0's
    """

    while len(seq) != max_len:
        seq.append(0)
    return seq

def loadMatrix(afile):
    alignment_matrix = []
    longest_len = 0
    matrix_file = open(afile)
    for line in matrix_file:
        alignment_matrix.append([])
        for item in line.strip().split():
            alignment_matrix[-1].append(int(item))
        seq_len = len(alignment_matrix[-1])
        if seq_len > longest_len:
            longest_len = seq_len    
    return [extendSeq(longest_len, seq) for seq in alignment_matrix]

def mkProteinPlot_2proteins(motif_name_file1, motif_name_file2,
                            motif_matrix_dir,
                            output_file):
    """ This makes the multiple alignment annotation figure for two proteins,
        one per column.

    @param motif_ls_file1: file with the locations and names of the motifs
                           separated by tabs, one per line
    @param motif_matrix_dir: directory with multiple 
                              alignments annotated w/ 0/1 for absence/presense
    @param motif_ls_file2: file with the locations and names of the motifs
                            separated by tabs, one per line   
    @param out_file: put the figure here
    """

    motifs1 = {}
    afile = open(motif_name_file1)
    for line in afile:
        [location, name] = line.strip().split('\t')
        motifs1[name] = location
        
    motifs2 = {}
    afile = open(motif_name_file2)
    for line in afile:
        [location, name] = line.strip().split('\t')
        motifs2[name] = location

    cols = 2
    rows = max([len(motifs1.keys()), len(motifs2.keys())])
    
    motif_index = 1
    #motif_figure = pylab.gcf()
    #default_size = motif_figure.get_size_inches()
    #motif_figure.set_size_inches( (default_size[0]*3, default_size[1]*3) )
    for motif in motifs1.keys():        
        motif_matrix = loadMatrix(motif_matrix_dir + motifs1[motif])
        #print rows, cols, motif_index
        pylab.subplot(rows, cols, motif_index)
        pylab.imshow(motif_matrix)
        pylab.xticks([], [])
        pylab.yticks([], [])
        pylab.title(motif)
        motif_index += 2
    motif_index = 2
    for motif in motifs2.keys():        
        motif_matrix = loadMatrix(motif_matrix_dir + motifs2[motif])
        #print rows, cols, motif_index
        pylab.subplot(rows, cols, motif_index)
        pylab.imshow(motif_matrix)
        pylab.xticks([], [])
        pylab.yticks([], [])
        pylab.title(motif)
        motif_index += 2
    pylab.savefig(output_file, dpi=(300))
    pylab.close()            

def mkProteinPlot(motif_name_file, motif_matrix_dir, output_file):
    """ This makes the multiple alignment annotation figure.

    @param motif_ls_file: file with the locations and names of the motifs
                          separated by tabs, one per line
    @param motif_matrix_dir: directory with multiple 
                             alignments annotated w/ 0/1 for absence/presense
    @param out_file: put the figure here
    """
    
    motifs = {}
    afile = open(motif_name_file)
    for line in afile:
        [location, name] = line.strip().split('\t')
        motifs[name] = location
    motif_count = len(motifs.keys())
    #print motifs
    if motif_count < 4:
        # up to 3x1
        cols = 1
        rows = motif_count
    elif motif_count < 7:
        # up to 3x2
        cols = 2
        rows = int(math.ceil( float(motif_count) / float(2) ))
    elif motif_count < 10:
        # up to 3x3
        cols = 3
        rows = int(math.ceil( float(motif_count) / float(3) ))
    elif motif_count < 13:
        # up to 3x4
        cols = 4
        rows = int(math.ceil( float(motif_count) / float(4) ))
    elif motif_count < 17:
        # up to 4x4
        rows = 4
        cols = int(math.ceil( float(motif_count) / float(4) ))
    elif motif_count < 21:
        # up to 5x4
        cols = 4
        rows = int(math.ceil( float(motif_count) / float(4) ))
    #cols = math.ceil( float(len(motifs.keys()))/float(5) )
    motif_index = 1
    #motif_figure = pylab.gcf()
    #default_size = motif_figure.get_size_inches()
    #motif_figure.set_size_inches( (default_size[0]*3, default_size[1]*3) )
    for motif in motifs.keys():        
        motif_matrix = loadMatrix(motif_matrix_dir + motifs[motif])
        #print rows, cols, motif_index
        pylab.subplot(rows, cols, motif_index)
        pylab.imshow(motif_matrix)
        pylab.xticks([], [])
        pylab.yticks([], [])
        pylab.title(motif)
        motif_index += 1
    pylab.savefig(output_file, dpi=(300))
    pylab.close()
