import MySQLdb, sys

# define species references in GO database
speciesID = dict()
speciesID['yeast'] = "97048"
speciesID['fly'] = "157890"
speciesID['human'] = '92001'
gocats = ['Process', 'Component', 'Function']

gocatAlias = dict()
gocatAlias['cellular_component'] = 'Component'
gocatAlias['molecular_function'] = 'Function'
gocatAlias['biological_process'] = 'Process'
gocatAlias['universal'] = 'Universal'

unknownGOterms = dict()
for t in ['GO:0005554', 'GO:0008372', 'GO:0000004']:
    unknownGOterms[t] = True

# connect to the GO database
def getGOConnection_Cursor():
    connection = MySQLdb.connect(host="mysql.ebi.ac.uk", \
                                 user='go_select', \
                                 db='go_latest', \
                                 port=4085, \
                                 passwd='amigo')
    cursor = connection.cursor()
    return [connection, cursor]

# for a given species,
#  for each GO term for that species
#   return a dict of the term's proteins and
#          an empty dict for the term's children
#          an empty dict for the term's parents
def getTermProperties(species):
    ls = dict()
    try:
        [connection, cursor] = getGOConnection_Cursor()
        cursor.execute("SELECT term_id, gene_product_id from association a inner join gene_product b on a.gene_product_id=b.id where b.species_id=" + speciesID[species] + ";")
        data = cursor.fetchall()
        for item in data:
            if not ls.has_key(item[0]):
                ls[ item[0] ] = [dict(), dict(), dict()] #proteins, children
            ls[ item[0] ][0][ item[1] ] = True
    except MySQLdb.OperationalError, message:
        print message[0]
        print message[1]
        sys.exit()

    cursor.close()
    connection.close()
    return ls

# get the parent/child relations
# for all terms
# return 2 lists
#  child2Parent and parent2Child
def getTermRelations():
    try:
        [connection, cursor] = getGOConnection_Cursor()
        # in term table
        #term2 (pair[1]) is child
        #term1 (pair[0]) is parent
        cursor.execute("select term1_id, term2_id, relationship_type_id from term2term;")
        data = cursor.fetchall()
        child2Parent = dict()
        parent2Child = dict()
        for pair in data:
            if not child2Parent.has_key(pair[1]):
                child2Parent[ pair[1] ] = dict()
            child2Parent[pair[1]][ pair[0] ] = True

            if not parent2Child.has_key(pair[0]):
                parent2Child[ pair[0] ] = dict()
            parent2Child[ pair[0] ][ pair[1] ] = True
            
    except MySQLdb.OperationalError, message:
        print message[0]
        print message[1]
        sys.exit()

    cursor.close()
    connection.close()
    return [child2Parent, parent2Child]

# return a dictionary and a maximum distance
#  from the root (dis from root to root is 0)
# indexed by distances from the root terms
# values are the GO terms found at a
#  specific distance
# a term can be at multiptle distances
#  from the root b/c there can be multiple paths
#  to a term
def getDistanceFromRoot():
    try:
        [connection, cursor] = getGOConnection_Cursor()
        cursor.execute("SELECT term.id, p.distance FROM term INNER JOIN graph_path AS p ON (p.term2_id=term.id) INNER JOIN term AS root ON (p.term1_id=root.id) WHERE root.is_root=1;")
        data = cursor.fetchall()
        disDict = dict()
        maxDist = 0
        for item in data:
            if not disDict.has_key(str(item[1])):
                if int(item[1]) > maxDist:
                    maxDist = int(item[1])
                disDict[str(item[1])] = dict()
            disDict[str(item[1])][item[0]] = True
    except MySQLdb.OperationalError, message:
        print message[0]
        print message[1]
        sys.exit()
        
    cursor.close()
    connection.close()
    return [disDict, maxDist]

# assigns children to all terms in termProperties
#  and assigns proteins and children to terms
#  that can be reached using child2Parent and a term
#  in termProperties
# GO doesn't give proteins for all its terms
#  so I must find this information using the
#  hierarchy
def fillInTermProperties(termProperties, child2Parent):
    distInfo = getDistanceFromRoot()
    distDict = distInfo[0]
    maxDist = distInfo[1]
    for dis in xrange(maxDist+1):
        for term in distDict[str(maxDist-dis)].keys():
            if termProperties.has_key(term):
                proteins = termProperties[term][0].keys()
                if child2Parent.has_key(term):
                    for p in child2Parent[term].keys():
                        if not termProperties.has_key(p):
                            termProperties[p] = [dict(), dict(), dict()]
                        for k in proteins:
                            termProperties[p][0][k] = True
                        termProperties[p][1][term] = True
                        for c in termProperties[term][1].keys():
                            termProperties[p][1][c] = True
            else:
                termProperties[term] = [dict(),dict(), dict()]
                if child2Parent.has_key(term):
                    for p in child2Parent[term].keys():
                        if not termProperties.has_key(p):
                            termProperties[p] = [dict(), dict(), dict()]
                        termProperties[p][1][term] = True    

# assigns children and parents to all terms in termProperties
#  and assigns proteins and children to terms
#  that can be reached using child2Parent or parent2Child
#  and a term in termProperties
# GO doesn't give proteins for all its terms
#  so I must find this information using the
#  hierarchy
def fillInTermPropertiesWParents(termProperties, child2Parent, parent2Child):
    fillInTermProperties(termProperties, child2Parent)
    for term in termProperties.keys():
        children = termProperties[term][1]
        for c in children.keys():
            termProperties[c][2][term] = True
                
# for each GO term id
# get its GOterm and cocat
def getTermIDtoGOterm():
    try:
        [connection, cursor] = getGOConnection_Cursor()
        cursor.execute(" SELECT id, acc, term_type, name FROM term WHERE term_type='biological_process' or term_type='cellular_component' or term_type='molecular_function' or term_type='universal';")
        data = cursor.fetchall()        
        termInfo = dict()
        for item in data:
            termInfo[ item[0] ] = [item[1], gocatAlias[item[2]], item[3]]
    except MySQLdb.OperationalError, message:
        print message[0]
        print message[1]
        sys.exit()
    cursor.close()
    connection.close()
    return termInfo        

# given 2 dicts of terms,
#  find the intersection,
#  split the intersecion into 3 gocats and print
def mkGOStandardFile(terms1, terms2, afile):
    ls = dict()
    for y in terms1.keys():
        if terms2.has_key(y):
            ls[y]=True
    termInfo = getTermIDtoGOterm()
    files = dict()
    for gocat in gocats:
        files[gocat] = open(afile+'.'+gocat, 'w')
    for item in ls.keys():
        itemName = termInfo[item][0]
        itemCat = termInfo[item][1]
        itemDesc = termInfo[item][2]
        files[itemCat].write(itemName + '\t' + itemDesc + '\n')
    for gocat in gocats:
        files[gocat].close()

def getTerm_proteins_children_parents(species):
    [child2Parent, parent2Child] = getTermRelations()
    termProperties = getTermProperties(species)
    fillInTermPropertiesWParents(termProperties, child2Parent, parent2Child)
    return termProperties

def getTermDesc(term):
    try:
        [connection, cursor] = getGOConnection_Cursor()
        cursor.execute("SELECT name from term WHERE id='" + \
                       term + "';")
        desc = cursor.fetchall()[0][0]
        return desc
    except:
        print 'error in getTermDesc: term =', term
        sys.exit(0)

def getTermSynonym(term):
    try:
        [connection, cursor] = getGOConnection_Cursor()
        cursor.execute("SELECT acc FROM term INNER JOIN term_synonym ON (term.id=term_synonym.term_id) WHERE acc_synonym='" + \
                       term + "';")
        data = cursor.fetchall()
        ls = []
        for item in data:
            ls.append(item[0])
        return ls
    except:
        print 'error in getTermSynonym: term =', term
        sys.exit(0)
        
