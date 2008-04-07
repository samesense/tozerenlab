from SOAPpy import WSDL


def printEdges(species_code):
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    pathways = serv.list_pathways(species_code)
    for path in pathways:
        path = path['entry_id']
        elements = serv.get_elements_by_pathway(path)
        relations = serv.get_element_relations_by_pathway(path)
        elementid2gene = dict()
        for element in elements:
            if element['type'] == 'gene':
                elementid = element['element_id']
                names = element['names']
                elementid2gene[elementid] = dict()
                for n in names:
                    elementid2gene[elementid][n] = True

        for relation in relations:
            k1 = relation['element_id1']
            k2 = relation['element_id2']
            if elementid2gene.has_key(k1) and elementid2gene.has_key(k2): 
                elements1 = elementid2gene[ k1 ]
                elements2 = elementid2gene[ k2 ]
                for e1 in elements1.keys():
                    for e2 in elements2.keys():
                        for s in relation['subtypes']:
                            print e1.split(':')[1] + '\t' + \
                                  e2.split(':')[1] + '\t' + \
                                  relation['type'] + '\t' + \
                                  s['relation'] + '\t' + path

def printPathwayDesc(species_code):
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    pathways = serv.list_pathways(species_code)
    for path in pathways:
        print path['entry_id'] + '\t' + path['definition']

def printGenes2Path(species_code):
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    pathways = serv.list_pathways(species_code)
    for path in pathways:
        genes = serv.get_genes_by_pathway(path[entry_id])
        for g in genes:
            print g + '\t' + path[entry_id]

def getPathwayEdges(path):
    edges = dict()
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    elements = serv.get_elements_by_pathway(path)
    relations = serv.get_element_relations_by_pathway(path)
    elementid2gene = dict()
    for element in elements:
        if element['type'] == 'gene':
            elementid = element['element_id']
            names = element['names']
            elementid2gene[elementid] = dict()
            for n in names:
                elementid2gene[elementid][n] = True

    for relation in relations:
        k1 = relation['element_id1']
        k2 = relation['element_id2']
        if elementid2gene.has_key(k1) and elementid2gene.has_key(k2): 
            elements1 = elementid2gene[ k1 ]
            elements2 = elementid2gene[ k2 ]
            for e1 in elements1.keys():
                for e2 in elements2.keys():
                    for s in relation['subtypes']:
                        gene1 = e1.split(':')[1]
                        gene2 = e2.split(':')[1]
                        if not edges.has_key(gene1):
                            edges[gene1] = dict()
                        if not edges.has_key(gene2):
                            edges[gene2] = dict()
                        if not edges[gene1].has_key(gene2):
                            edges[gene1][gene2] = dict()
                        if not edges[gene2].has_key(gene1):
                            edges[gene2][gene1] = dict()
                        edges[gene1][gene2][ relation['type'] + '/' + s['relation'] ] = True
                        edges[gene2][gene1][ relation['type'] + '/' + s['relation'] ] = True
    return edges

def getPathwayGenes(path):
    genes = dict()
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    elements = serv.get_elements_by_pathway(path)
    for element in elements:
        if element['type'] == 'gene':            
            names = element['names']            
            for n in names:
                genes[n.split(':')[1]] = True
    return genes

def getPathwayComplexDict(path):
    complex2gene = dict()
    gene2complex = dict()
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    elements = serv.get_elements_by_pathway(path)
    for element in elements:
        if element['type'] == 'gene':            
            names = element['names']
            elementID = str(element['element_id'])
            complex2gene[elementID] = dict()
            for n in names:
                n=n.split(':')[1]
                complex2gene[elementID][n] = True
                if not gene2complex.has_key(n): gene2complex[n] = dict()
                gene2complex[n][elementID] = True
    return [complex2gene, gene2complex]

