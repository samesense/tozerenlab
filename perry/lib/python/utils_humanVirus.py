#---------------------------------------
#
# Author     Perry Evans
#            evansjp@mail.med.upenn.edu
# 2008
#---------------------------------------
"""
This holds various helper functions for the 
Human_Virus project.
"""
import string, sys
sys.path.append('/home/perry/Code/lib/python/')
import utils_graph, utils_motif

def initResults(results_dict):
    """ Initialize this ELM entry with 
        empty results {}.

    @param results_dict: {} to init
    """

    results_dict['H1'] = {}
    results_dict['H2'] = {}
    results_dict['H12'] = {}
    results_dict['H1toH2'] = {}

def recordResults(results, gene1, gene2):
    """ Fill in the results for a motif search
        for an elm.

    @param results: {} with h1, h2, h12, h1toh2
                    entries.
    @param gene1: protein id for node.
    @param gene2: protein id for node.
    """

    results['H1'][gene1] = True
    results['H2'][gene2] = True
    results['H12'][gene1] = True
    results['H12'][gene2] = True
    if not results['H1toH2'].has_key(gene1):
        results['H1toH2'][gene1] = {}
    results['H1toH2'][gene1][gene2] = True

def recordMethodResults(viral_protein, elm, gene1, gene2,
                        viral_protein_results, viralELMs,
                        all_virus_results, elm_results, dir):
    """ For a method like Network, fill in all the {}
        I'm tracking to update this predicted edge.

    @param viral_protein: 
    @param elm: 
    @param gene1:
    @param gene2:
    @param viral_protein_results:
    @param viralELMs:
    @param all_virus_resutls:
    @param elm_results:
    @param dir
    """

    if not elm_results.has_key(elm):
        elm_results[elm] = {}
        initResults(elm_results[elm])
    if not viralELMs.has_key(elm):
        viralELMs[elm] = {}
        initResults(viralELMs[elm])
    recordResults(viralELMs[elm], gene1, gene2)
    recordResults(elm_results[elm], gene1, gene2)         
    recordResults(viral_protein_results, gene1, gene2)
    recordResults(all_virus_results, gene1, gene2)
    all_virus_results['virus_H1_H2'][viral_protein + '>' + gene1 + '>'
                                      + gene2] = True
def dumpResults(results, dir, prefix):
    """ Print out the results {} in 4 files,
        using this prefix & directory.

    @parm results: {} w/ H1, H2, H12 {} and edge {} H1toH2.
    @param dir: directory to put results, /w '/' at end.
    @param prefix: name of file to stick before ex. '.H1'.
    """

    utils_graph.dumpNodes(dir + prefix + '.H1', results['H1'])
    utils_graph.dumpNodes(dir + prefix + '.H2', results['H2'])
    utils_graph.dumpNodes(dir + prefix + '.H12', results['H12'])
    utils_graph.dumpEdges(dir + prefix + '.H1toH2', results['H1toH2'])

def dumpMotifSearchResults(viral_protein, elm, viralELMs,
                           viral_protein_results, dir):
    dumpResults(viral_protein_results, dir, viral_protein)
    for elm in viralELMs.keys():
        dumpResults(viralELMs[elm], dir, viral_protein + '.' + elm)

def wrapMotifSearchDump(all_virus_results, dir):
    """
    
    @param all_virus_results:
    @param dir: 
    """

    dumpResults(all_virus_results, dir, 'all')
    virus_edge_file = open(dir + 'virus_h1_h2', 'w')
    for entry in all_virus_results['virus_H1_H2'].keys():
        [viral_protein, gene1, gene2] = entry.split('>')
        virus_edge_file.write(viral_protein + '\t'
                              + gene1 + '\t' + gene2 + '\n')   
    virus_edge_file.close()

def cp2(elm_ls, new_ls):
    for elm in elm_ls.keys():
        if not new_ls.has_key(elm): new_ls[elm] = {}
        for g in elm_ls[elm].keys():
            new_ls[elm][g] = True

def union(elmLs1, elmLs2):
    newLs = {}
    cp2(elmLs1, newLs)
    cp2(elmLs2, newLs)
    return newLs

def get_virus2elm(conservedELMsFile):
    virus2elm = {}
    f = open(conservedELMsFile)
    for line in f.xreadlines():
        [v, elm] = map(string.strip, line.split('\t'))
        if not virus2elm.has_key(v):
            virus2elm[v] = {}
        virus2elm[v][elm] = {}
    f.close()
    return virus2elm

def getViralELMs(conservedELMsFile):
    elms = {}
    f = open(conservedELMsFile)
    for line in f.xreadlines():
        [v, elm] = map(string.strip, line.split('\t'))
        elms[elm] = {}
    f.close()
    return elms

def getProteinsForELMs(elm2protein_file, useELMs, proteins):
    f = open(elm2protein_file)
    for line in f.xreadlines():
        if line[0] != '#':
            [elm, gene] = map(string.strip, line.split())
            if useELMs.has_key(elm):
                useELMs[elm][gene] = True
                if not proteins.has_key(gene):
                    proteins[gene] = {}
                proteins[gene][elm] = True
    f.close()

def expandProteinsForELMs(humanAnnotationFile, elm_pairs_dir, 
                          useELMs, proteins, domain_tools):
    prosite2protein = utils_motif.annotation2protein(humanAnnotationFile,
                                                     domain_tools)
    for domain in domain_tools.keys():
        f = open(elm_pairs_dir + 'ELM.' + domain + '.pairs')
        for line in f.readlines():
            sp = map(string.strip, line.split('\t'))
            if len(sp) > 1:
                if sp[1] != '':
                    elm = sp[0]
                    if useELMs.has_key(elm):
                        for protein in prosite2protein[ sp[1] ]:
                            useELMs[elm][protein] = True
                            if not proteins.has_key(protein): proteins[protein] = {}
                            proteins[protein][elm] = True                 
        f.close()

def predict(useELMs, domain_tools, netFile):
    h1 = {}
    h2 = {}
    net = utils_graph.getEdges(netFile)
    proteins = {}
    tool_d = {}
    tool_d['ELM'] = True
    protein2elm = utils_motif.protein2annotation('/home/perry/Projects/Human_Virus/Data/human.annotations', tool_d)
    getProteinsForELMs(useELMs, proteins)
    expandProteinsForELMs(useELMs, proteins, domain_tools)
    for g1 in proteins.keys():
        if net.has_key(g1):
            for elm in proteins[g1].keys():
                for g2 in protein2elm.keys():
                    if net[g1].has_key(g2) and protein2elm[g2].has_key(elm) and g1 != g2:
                        if not h1.has_key(elm): h1[elm] = {}
                        if not h2.has_key(elm): h2[elm] = {}
                        h1[elm][g1] = True
                        h2[elm][g2] = True
    return [h1, h2]

def convert2dict(elm_ls):
    elm_d = {}
    for elm in elm_ls.keys():
        elm_d[elm] = {}
    return elm_d

def H1(useELMs, domain_tools, netFile):
    elm_d = convert2dict(useELMs)
    return predict(elm_d, domain_tools, netFile)[0]
    
def H2(useELMs, domain_tools, netFile):
    elm_d = convert2dict(useELMs)
    return predict(elm_d, domain_tools, netFile)[1]

def H12(useELMs, domain_tools, netFile):
    elm_d = convert2dict(useELMs)
    [h1, h2] = predict(elm_d, domain_tools, netFile)
    return intersect(h1, h2)

def predictions(useELMs, domain_tools, netFile):
    elm_d = convert2dict(useELMs)
    [h1, h2] = predict(elm_d, domain_tools, netFile) 
    return [h1, h2, union(h1, h2)] 
