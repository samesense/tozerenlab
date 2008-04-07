import string

def getNodes(afile):
    nodes = dict()
    f = open(afile)
    for line in f.xreadlines():
        nodes[ line.strip() ] = True
    f.close()
    return nodes

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
