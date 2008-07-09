import string

def protein2motif(afile):
    d = dict()
    f = open(afile)
    for line in f.xreadlines():
        sp = line.split('\t')
        if not d.has_key( sp[0] ): d[ sp[0] ] = {}
        d[ sp[0] ][ sp[3] ] = True
    f.close()
    return d

def motif2protein(afile):
    d = dict()
    f = open(afile)
    for line in f.xreadlines():
        sp = line.split('\t')
        if not d.has_key( sp[3] ): d[ sp[3] ] = {}
        d[ sp[3] ][ sp[0] ] = True
    f.close()
    return d

def annotation2protein(afile, tool_dict):
    d = dict()
    f = open(afile)
    for line in f.xreadlines():
        [geneid, start, stop,
         domain, seq_desc, tool] = map(string.strip, line.split('\t'))
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
    d = dict()
    f = open(afile)
    for line in f.xreadlines():
        [geneid, start, stop,
         domain, seq_desc, tool] = map(string.strip, line.split('\t'))
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
         domain, seq_desc, tool] = map(string.strip, line.split('\t'))
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
         domain, seq_desc, tool] = map(string.strip, line.split('\t'))
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
        
def addBackslash(s):
    return s.replace('_', '\_')
