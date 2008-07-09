#----------------------------------------
#
# Author     Perry Evans
#            evansjp@mail.med.upenn.edu
# 2008
#----------------------------------------
"""
Functions here help with network loading,
printing, and searching.
"""
import string

def getNodes(afile):
    """ Make a {} from a file.  Each line
        of the file is one item in the {}.
    
    @param afile: file w/ one item per line
    @return: {} of nodes
    """

    nodes = dict()
    f = open(afile)
    for line in f:
        nodes[ line.strip() ] = True
    f.close()
    return nodes

def dumpNodes(afile, ls):
    """ Output {} to file.  One item
    per line.

    @param afile: write to this file
    @param ls: {} with items for writing
    """

    f = open(afile, 'w')
    for g in ls.keys():
        f.write(g + '\n')
    f.close()

def getEdges(afile):
    """ Make edge {} from file.
    {} format is d[n1][n2] = True

    @param afile: tab-delimited file; one edge per line; edges are bi-directional
    @return: {} of edges
    """

    edges = dict()
    f = open(afile)
    for line in f():
        [n1, n2] = map(string.strip, line.split('\t'))
        if not edges.has_key(n1): edges[n1] = dict()
        if not edges.has_key(n2): edges[n2] = dict()
        edges[n1][n2] = True
        edges[n2][n1] = True
    f.close()
    return edges

def dumpEdges(afile, d):
    """ Write edges to file.
    
    @param afile: output is tab-delimited.
    bi-directional edges are not duplicated
    @param d: {} format is d[n1][n2] = True
    """
    
    f = open(afile, 'w')
    seen = {}
    for g1 in d.keys():
        for g2 in d[g1].keys():
            k1 = g1+':'+g2
            k2 = g2+':'+g1
            if not seen.has_key(k1) and not seen.has_key(k2):
                seen[k1] = True
                seen[k2] = True
                f.write(g1 + '\t' + g2 + '\n')
    f.close()

def mkNodesFromEdges(afile):
    """ Given an edge file, make a {} of nodes.

    @param afile: tab-delimted file of edges; one per line; edges are bi-directional
    @return: {} of nodes
    """

    nodes = {}
    f = open(afile)
    for line in f.xreadlines():
        sp = line.split('\t')
        nodes[sp[0]] = True
        nodes[sp[1]] = True
    f.close()
    return nodes

def getConnected(ls_file, net_file):
    """ Find all the genes in the network
    that are one away from these genes.

    @param ls_file: file of genes to find connections to; one gene per line
    @param net_file: tab-delimited file for the network; one edge per line; edges are bi-directional
    @return: {} of genes
    """

    ls = getNodes(ls_file)
    net = getEdges(net_file)
    second_level_ls = {}
    for gene in ls.keys():
        if net.has_key(gene):
            for g2 in net[gene].keys():
                if gene != g2:
                    second_level_ls[g2] = True
    return second_level_ls

def intersectLists(list_of_lists):
    """ Given a list of {}, return the intersection
    as a {}.

    @param list_of_lists: [] of {}
    @return: {} of genes
    """

    same = {}
    for s in list_of_lists[0].keys():
        count = 0
        for ls in list_of_lists[1:]:
            if ls.has_key(s): count += 1
        if count == len(list_of_lists)-1:
            same[s] = True
    return same

def unionLists(list_of_lists):
    """ Given a list of {}, return the union
    as a {}.

    @param list_of_lists: list of {}
    @return: {} of union
    """

    all = {}
    for ls in list_of_lists:
        for k in ls.keys():
            all[k] = True
    return all

def degree(gene, network):
    """ Find the degree of this node
    in this network.
    
    @param gene: node in the network
    @param network: edge {} format d[n1][n2] = True
    @return: deg of gene as int
    """

    return len(network[gene].keys())

def avgConnectivity(geneLs, network):
    """ Given a set of genes, return
    their average connectivity in the network.

    @param geneLs: {} of genes
    @param network: edge {} format d[n1][n2] = True
    @return: float for avg connectivity
    """

    [connectivity, total] = [0, 0]
    for gene in geneLs.keys():
        if network.has_key(gene):
            connectivity += degree(gene, network)
            total += 1
    return float(connectivity) / float(total)
