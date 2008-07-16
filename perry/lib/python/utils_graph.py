#----------------------------------------
#
# Author     Perry Evans
#            evansjp@mail.med.upenn.edu
# 2008
#
#----------------------------------------
"""
Functions here help with network loading,
printing, and searching.
"""

def getNodes(afile):
    """ Make a {} from a file.  Each line
        of the file is one item in the {}.
    
    @param afile: file w/ one item per line
    @return: {} of nodes
    """

    nodes = dict()
    f_nodes = open(afile)
    for line in f_nodes:
        nodes[ line.strip() ] = True
    f_nodes.close()
    return nodes

def dumpNodes(afile, node_dict):
    """ Output {} to file.  One item
    per line.

    @param afile: write to this file
    @param node_dict: {} with items for writing
    """

    f_out = open(afile, 'w')
    for node in node_dict.keys():
        f_out.write(node + '\n')
    f_out.close()

def getEdges(afile):
    """ Make edge {} from file.
    {} format is d[n1][n2] = True

    @param afile: tab-delimited file; one edge per line; edges are bi-directional
    @return: {} of edges
    """

    edges = dict()
    f_edges = open(afile)
    for line in f_edges:
        [node1, node2] = [x.strip() for x in line.split('\t')]
        if not edges.has_key(node1): 
            edges[node1] = dict()
        if not edges.has_key(node2): 
            edges[node2] = dict()
        edges[node1][node2] = True
        edges[node2][node1] = True
    f_edges.close()
    return edges

def dumpEdges(afile, edge_dict):
    """ Write edges to file.
    
    @param afile: output is tab-delimited; bi-directional edges are not duplicated
    @param edge_dict: {} format is d[n1][n2] = True
    """
    
    f_out = open(afile, 'w')
    seen = {}
    for gene1 in edge_dict.keys():
        for gene2 in edge_dict[gene1].keys():
            key1 = gene1+':'+gene2
            key2 = gene2+':'+gene1
            if not seen.has_key(key1) and not seen.has_key(key2):
                seen[key1] = True
                seen[key2] = True
                f_out.write(gene1 + '\t' + gene2 + '\n')
    f_out.close()

def mkNodesFromEdges(afile):
    """ Given an edge file, make a {} of nodes.

    @param afile: tab-delimted file of edges; one per line; edges are bi-directional
    @return: {} of nodes
    """

    nodes = {}
    f_edges = open(afile)
    for line in f_edges:
        splitup = line.split('\t')
        nodes[splitup[0]] = True
        nodes[splitup[1]] = True
    f_edges.close()
    return nodes

def getConnected(ls_file, net_file):
    """ Find all the genes in the network
    that are one away from these genes.

    @param ls_file: file of genes to find connections to; one gene per line
    @param net_file: tab-delimited file for the network; one edge per line; edges are bi-directional
    @return: {} of genes
    """

    seed_ls = getNodes(ls_file)
    net = getEdges(net_file)
    second_level_ls = {}
    for gene in seed_ls.keys():
        if net.has_key(gene):
            for gene2 in net[gene].keys():
                if gene != gene2:
                    second_level_ls[gene2] = True
    return second_level_ls

def intersectLists(list_of_lists):
    """ Given a list of {}, return the intersection
    as a {}.

    @param list_of_lists: [] of {}
    @return: {} of genes
    """

    same = {}
    for item in list_of_lists[0].keys():
        count = 0
        for alist in list_of_lists[1:]:
            if alist.has_key(item): 
                count += 1
        if count == len(list_of_lists)-1:
            same[item] = True
    return same

def intersectFiles(list_of_files):
    """ Given a list of files, return the intersection
    as a {}.

    @param list_of_lists: [] of files
    @return: {} of genes
    """

    list_of_lists = []
    for afile in list_of_files:
        list_of_lists.append(getNodes(afile))
    return intersectLists(list_of_lists)
    
def unionLists(list_of_lists):
    """ Given a list of {}, return the union
    as a {}.

    @param list_of_lists: list of {}
    @return: {} of union
    """

    union = {}
    for alist in list_of_lists:
        for k in alist.keys():
            union[k] = True
    return union

def degree(gene, network):
    """ Find the degree of this node
    in this network.
    
    @param gene: node in the network
    @param network: edge {} format d[n1][n2] = True
    @return: deg of gene as int
    """

    return len(network[gene].keys())

def avgConnectivity(gene_ls, network):
    """ Given a set of genes, return
    their average connectivity in the network.

    @param geneLs: {} of genes
    @param network: edge {} format d[n1][n2] = True
    @return: float for avg connectivity
    """

    [connectivity, total] = [0, 0]
    for gene in gene_ls.keys():
        if network.has_key(gene):
            connectivity += degree(gene, network)
            total += 1
    return float(connectivity) / float(total)

def getDegreeDistribution(gene_dict, network_edges):
    """ For these genes, return the degree distribution
        according to this network.

    @param gene_dict: {} w/ gene names
    @param network_edges {} of edges d[n1][n2] = True
    @return {} of degrees d[degree] = count
    """

    degreeDict = []
    for gene in gene_dict.keys():
        if network_edges.has_key(gene):
            degreeDict.append(degree(gene, network_edges))
            #if not degreeDict.has_key(deg):
            #    degreeDict[str(deg)] = 0
            #degreeDict[str(deg)] += 1
    return degreeDict
