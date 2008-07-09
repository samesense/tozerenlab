import sys#, urllib2, ClientForm, httplib, ClientCookie, time
sys.path.append('/home/perry/Code/lib/python/will/')
import PyMozilla
import string, re
#import utils_kegg

def fullMonty(s):
    if s == '': return s
    if s[-1] == 'O':
        s = s[0:-4]
    elif s[-8:-5] == 'hsa':
        s = s[0:-10]
    return s.strip().lstrip().strip(',')

def cat2term2protein_sigFile(afile):
    d = {}
    f = open(afile)
    f.readline()
    for line in f.xreadlines():
        [Category, Term, Count, Percent,
         PValue, Genes, ListTotal,
         PopHits, PopTotal, FoldChange] = map(string.strip, line.split('\t'))
        if not d.has_key(Category):
            d[Category] = {}
        genes = Genes.strip(',').split(',')
        geneLs = {}
        for gene in genes: geneLs[gene.strip().lstrip()] = True
        d[Category][Term] = {}
        d[Category][Term]['count'] = int(Count)
        d[Category][Term]['pvalue'] = float(PValue)
        d[Category][Term]['genes'] = geneLs
    f.close()
    return d

def cat2term2protein_sigLs(ls):
    d = {}
    for line in ls[1:-1]:
        [Category, Term, Count, Percent,
         PValue, Genes, ListTotal,
         PopHits, PopTotal, FoldChange] = map(string.strip, line.split('\t'))
        if not d.has_key(Category):
            d[Category] = {}
        genes = Genes.strip(',').split(',')
        geneLs = {}
        for gene in genes: geneLs[gene.strip().lstrip()] = True
        d[Category][Term] = {}
        d[Category][Term]['count'] = int(Count)
        d[Category][Term]['pvalue'] = float(PValue)
        d[Category][Term]['genes'] = geneLs
    return d

def cat2term2protein_allFile(afile):
    d = {}
    f = open(afile)
    Categories = map(string.strip, f.readline().split('\t')[3:7])
    for cat in Categories: d[cat] = {}
    for line in f.xreadlines():
        sp = map(string.strip, line.split('\t'))
        for i in xrange(len(Categories)):
            cat = Categories[i]
            terms = map(string.strip, sp[i+3].strip(',').split(':'))
            for x in xrange(len(terms)-1):
                if terms[x] != '':
                    term = terms[x].split()[-1]+':'+fullMonty(terms[x+1])
                    #print cat, term
                    if not d[cat].has_key(term):
                        d[cat][term] = {}
                        d[cat][term]['count'] = 0
                        d[cat][term]['genes'] = {}
                    d[cat][term]['genes'][sp[0]] = True
                    d[cat][term]['count'] += 1               
    f.close()
    return d

def cat2term2protein_allLs(ls):
    d = {}
    Categories = map(string.strip, ls[0].split('\t')[3:7])
    for cat in Categories: d[cat] = {}
    for line in ls[1:-1]:
        sp = map(string.strip, line.split('\t'))
        for i in xrange(len(Categories)):
            cat = Categories[i]
            terms = map(string.strip, sp[i+3].strip(',').split(':'))
            for x in xrange(len(terms)-1):
                if terms[x] != '':
                    term = terms[x].split()[-1]+':'+fullMonty(terms[x+1])
                    #print cat, term
                    if not d[cat].has_key(term):
                        d[cat][term] = {}
                        d[cat][term]['count'] = 0
                        d[cat][term]['genes'] = {}
                    d[cat][term]['genes'][sp[0]] = True
                    d[cat][term]['count'] += 1               
    return d

def getSigKEGGgenes_sigFile(afile, cutoff):
    genes = {}
    cat2term2protein = cat2term2protein_sigFile(afile)
    for path in cat2term2protein['KEGG_PATHWAY'].keys():
        if cat2term2protein['KEGG_PATHWAY'][path]['pvalue'] < cutoff:
            for g in cat2term2protein['KEGG_PATHWAY'][path]['genes'].keys():
                genes[g] = True
    return genes

def getSigKEGGgenes(motifSearchResults, bg, cutoff):
    
    genes = {}
    cat2term2protein = getSigTerms(motifSearchResults, bg)
    for path in cat2term2protein['KEGG_PATHWAY'].keys():
        if cat2term2protein['KEGG_PATHWAY'][path]['pvalue'] < cutoff:
            for g in cat2term2protein['KEGG_PATHWAY'][path]['genes'].keys():
                genes[g] = True
    return genes

#def getGenesForEnrichedKEGG(afile, cutOff):
#    genes = {}
#    cat2term2protein = cat2term2protein_sigFile(afile)
#    paths = {}
#    for kegg in cat2term2protein['KEGG_PATHWAY'].keys():
#        paths['path:' + kegg.split(':')[0]] = True
#    for kegg in paths.keys():
#        path_genes = utils_kegg.getPathwayGenes(kegg)
#        for gene in path_genes.keys():
#            genes[gene] = True
#    return genes

def getEnrichedTerms(fg, bg):
    pass

def getAllTerms(geneFile):
    MozEmu = PyMozilla.MozillaEmulator(cacher=None,trycount=0)
    #try:
    files = [['fileBrowser', 'testData', '']]
    fields = []
    #fields.append(['removeIndex', '0'])
    MozEmu.download('http://david.abcc.ncifcrf.gov/tools.jsp')
    #for i in xrange(10):
    #    MozEmu.post_multipart('http://david.abcc.ncifcrf.gov/tools.jsp', fields, files)
        #MozEmu.download('http://david.abcc.ncifcrf.gov/tools.jsp')
    
    fdata = file(geneFile).read()
    fields = []
    fields.append(['idType', 'ENTREZ_GENE_ID'])
    fields.append(['Mode', 'file'])
    fields.append(['uploadType','list'])

    files = []
    files.append(['fileBrowser', 'testData', fdata])
    
    ans = MozEmu.post_multipart('http://david.abcc.ncifcrf.gov/tools.jsp', fields, files)
    
    david = MozEmu.download('http://david.abcc.ncifcrf.gov/template.css')
    david = MozEmu.download('http://david.abcc.ncifcrf.gov/scripts/sidebar.js')
    david = MozEmu.download('http://david.abcc.ncifcrf.gov/scripts/summary.js')
    userDownloadPage = MozEmu.download('http://david.abcc.ncifcrf.gov/annotationReport.jsp?annot=,GOTERM_BP_5,GOTERM_CC_5,GOTERM_MF_5,KEGG_PATHWAY&currentList=0')
    index = userDownloadPage.find('UserDownload')
    addr = userDownloadPage[index+13:].split('.')[0]
    return cat2term2protein_allLs(MozEmu.download('http://david.abcc.ncifcrf.gov/UserDownload/' + addr  + '.txt').split('\n'))

def getSigTerms(fgFile, bgFile):
    MozEmu = PyMozilla.MozillaEmulator(cacher=None,trycount=0)
    #try:
    files = [['fileBrowser', 'testData', '']]
    fields = []
    #fields.append(['removeIndex', '0'])
    MozEmu.download('http://david.abcc.ncifcrf.gov/tools.jsp')
    #for i in xrange(10):
    #    MozEmu.post_multipart('http://david.abcc.ncifcrf.gov/tools.jsp', fields, files)
        #MozEmu.download('http://david.abcc.ncifcrf.gov/tools.jsp')
    
    fdata = file(fgFile).read()
    fields = []
    fields.append(['idType', 'ENTREZ_GENE_ID'])
    fields.append(['Mode', 'file'])
    fields.append(['uploadType','list'])

    files = []
    files.append(['fileBrowser', 'testData', fdata])
    
    ans = MozEmu.post_multipart('http://david.abcc.ncifcrf.gov/tools.jsp', fields, files)
    
    david = MozEmu.download('http://david.abcc.ncifcrf.gov/template.css')
    david = MozEmu.download('http://david.abcc.ncifcrf.gov/scripts/sidebar.js')
    david = MozEmu.download('http://david.abcc.ncifcrf.gov/scripts/summary.js')

    # upload the background
    fdata = file(bgFile).read()
    fields = []
    fields.append(['idType', 'ENTREZ_GENE_ID'])
    fields.append(['Mode', 'file'])
    fields.append(['uploadType','population'])

    files = []
    files.append(['fileBrowser', 'testData', fdata])
    
    ans = MozEmu.post_multipart('http://david.abcc.ncifcrf.gov/tools.jsp', fields, files)
    
    david = MozEmu.download('http://david.abcc.ncifcrf.gov/template.css')
    david = MozEmu.download('http://david.abcc.ncifcrf.gov/scripts/sidebar.js')
    david = MozEmu.download('http://david.abcc.ncifcrf.gov/scripts/summary.js')
    #david = MozEmu.download('http://david.abcc.ncifcrf.gov/scripts/summary.jsp')

    userDownloadPage = MozEmu.download('http://david.abcc.ncifcrf.gov/chartReport.jsp?annot=,GOTERM_BP_5,GOTERM_CC_5,GOTERM_MF_5,KEGG_PATHWAY&currentList=0')
    index = userDownloadPage.find('UserDownload')
    addr = userDownloadPage[index+13:].split('.')[0]
    res = MozEmu.download('http://david.abcc.ncifcrf.gov/UserDownload/' + addr  + '.txt')
    return cat2term2protein_sigLs(res.split('\n'))
