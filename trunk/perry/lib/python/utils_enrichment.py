#-----------------------------------------
#  Author:     Perry Evans
#              evansjp@mail.med.upenn.edu
#  2008
#
#-----------------------------------------
"""
Functions for p-values of enrichment of motifs
by permutation tests.
"""

import random
import utils_graph, utils_motif

def scoreProteinLs(alist, annotations, motifs):
    """ Given a list of genes, and annotations
        for these genes, count the # of times
        each given motif is seen.  Motifs can
        occur multiple times on a gene.

    @param alist: {} of genes
    @param annotations: {} annotations[protein][motif] = [motif_instances]
    @param motifs: {} of motifs
    @return: {} [motif] = # of times seen
    """

    scores = {}
    for motif in motifs.keys():
        scores[motif] = 0
    for protein in alist.keys():
        if annotations.has_key(protein):
            for motif in motifs.keys():
                if annotations[protein].has_key(motif):
                    scores[motif] += len(annotations[protein][motif])
    return scores

def scoreBckgnd(size, alist, annotations, motifs):
    """ Make a background permutation of this size by
        randomly drawing from the given list.  Return
        motif counts for this permutation.    
    
    @param size: # of background draws for this permutation
    @param alist: background pool
    @param annotations: {}[protein][motif] = [motif_instance]
    @param motifs: motifs to get background scores for
    @return: {} [motif] = # of times seen in this permutation
    """

    bck_ls = {}
    lenls = len(alist) - 1
    keys = alist.keys()
    while len(bck_ls) < size:
        rand_index = random.randint(0, lenls)
        bck_ls[ keys[rand_index] ] = True
    return scoreProteinLs(bck_ls, annotations, motifs)        

def findEnrichedAnnotations(fgnd_file, bgnd_file, annotation_file, motifs):
    """  Check to see if the given motifs are enriched in the fg,
         as compared to the bg.  Return motif {}, each entry
         is [#of proteins with motif, |fg|, # of times bg
         permutation is as good as fg].  Do 1000 background
         permutations.

    @param fgnd_file: file of foreground genes
    @param bgnd_file: file of backgound genes
    @param annotation_file: one annotation per line; 
    @param motifs: {} of motifs you are checking for enrichment
    @return: [#of proteins with motif, |fg|, # of times bg is as good as fg] for each motif
    """

    fg_ls = utils_graph.getNodes(fgnd_file)
    bg_ls = utils_graph.getNodes(bgnd_file)
    annotations = utils_motif.protein2annotation_forMotifs(annotation_file,
                                                           motifs)
    scores = scoreProteinLs(fg_ls, annotations, motifs)
    better_counts = {}
    for motif in motifs.keys():
        better_counts[motif] = 0
    for i in xrange(1000):
        bck_scores = scoreBckgnd( len(fg_ls),  bg_ls, annotations, motifs)
        for motif in motifs.keys():
            if bck_scores[motif] >= scores[motif]:
                better_counts[motif] += 1
    for motif in motifs.keys():
        scores[motif] = [ scores[motif], len(fg_ls), better_counts[motif] ]
    return scores
    
