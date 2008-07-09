import sys, random
sys.path.append('/home/perry/Code/lib/python/')
import utils_graph, utils_motif

def scoreProteinLs(ls, annotations, motifs):
    scores = {}
    for m in motifs.keys():
        scores[m] = 0
    for protein in ls.keys():
        if annotations.has_key(protein):
            for m in motifs.keys():
                if annotations[protein].has_key(m):
                    scores[m] += len(annotations[protein][m])
    return scores

def scoreBckgnd(size, ls, annotations, motifs):
    bck_ls = {}
    ll = len(ls)-1
    keys = ls.keys()
    while len(bck_ls) < size:
        r = random.randint(0, ll)
        bck_ls[ keys[r] ] = True
    return scoreProteinLs(bck_ls, annotations, motifs)        

def findEnrichedAnnotations(fgnd_file, bgnd_file, annotation_file, motifs):
    fg_ls = utils_graph.getNodes(fgnd_file)
    bg_ls = utils_graph.getNodes(bgnd_file)
    annotations = utils_motif.protein2annotation_forMotifs(annotation_file, motifs)
    scores = scoreProteinLs(fg_ls, annotations, motifs)
    better_counts = {}
    for m in motifs.keys():
        better_counts[m] = 0
    for i in xrange(1000):
        bck_scores = scoreBckgnd( len(fg_ls),  bg_ls, annotations, motifs)
        for m in motifs.keys():
            if bck_scores[m] >= scores[m]:
                better_counts[m] += 1
    for m in motifs.keys():
        scores[m] = [ scores[m], len(fg_ls), better_counts[m] ]
    return scores
    
