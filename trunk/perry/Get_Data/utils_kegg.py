from SOAPpy import WSDL
import sys

wsdl = 'http://soap.genome.jp/KEGG.wsdl'
serv = WSDL.Proxy(wsdl)
pathways = serv.list_pathways('hsa')
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
                        print e1.split(':')[1] + '\t' + e2.split(':')[1] + '\t' + relation['type'] + '\t' + s['relation'] + '\t' + path
