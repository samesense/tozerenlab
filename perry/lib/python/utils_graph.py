import string

def getNodes(afile):
    nodes = dict()
    f = open(afile)
    for line in f.xreadlines():
        nodes[ line.strip() ] = True
    f.close()
    return nodes

def dumpNodes(afile, ls):
    f = open(afile, 'w')
    for g in ls.keys():
        f.write(g + '\n')
    f.close()

def getEdges(afile):
    edges = dict()
    f = open(afile)
    for line in f.xreadlines():
        [n1, n2] = map(string.strip, line.split('\t'))
        if not edges.has_key(n1): edges[n1] = dict()
        if not edges.has_key(n2): edges[n2] = dict()
        edges[n1][n2] = True
        edges[n2][n1] = True
    f.close()
    return edges

def dumpEdges(afile, d):
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
    nodes = {}
    f = open(afile)
    for line in f.xreadlines():
        sp = line.split('\t')
        nodes[sp[0]] = True
        nodes[sp[1]] = True
    f.close()
    return nodes

def getConnected(ls_file, net_file):
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
    same = {}
    for s in list_of_lists[0].keys():
        count = 0
        for ls in list_of_lists[1:]:
            if ls.has_key(s): count += 1
        if count == len(list_of_lists)-1:
            same[s] = True
    return same

def unionLists(list_of_lists):
    all = {}
    for ls in list_of_lists:
        for k in ls.keys():
            all[k] = True
    return all

# get degree
def degree(gene, network):
    return len(network[gene].keys())

# given a set of genes
# return their average connectivity
def avgConnectivity(geneLs, network):
    [connectivity, total] = [0, 0]
    for gene in geneLs.keys():
        if network.has_key(gene):
            connectivity += degree(gene, network)
            total += 1
    return float(connectivity) / float(total)
