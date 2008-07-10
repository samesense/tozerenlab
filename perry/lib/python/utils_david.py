#-----------------------------------------
#  Author:     Perry Evans
#              evansjp@mail.med.upenn.edu
#  2008
#
#-----------------------------------------
"""
Functions for automating calls to DAVID
"""
import sys
import PyMozilla
import re

def fullMonty(astring):
    """ Remove commas and whitespace from
        this term.  This only works for
        human proteins.

    @param astring: string to be formatted
    @return: formatted string
    """

    if astring == '': return astring
    if astring[-1] == 'O':
        astring = astring[0:-4]
    elif astring[-8:-5] == 'hsa':
        astring = astring[0:-10]
    return astring.strip().lstrip().strip(',')

def cat2term2protein_sig(iterable):
    """ Take the output from DAVID chart
        and organize by category, term, protein.
    
    @param iterable: open file or list that has DAVID output
    @return: [Category][Term]['count' | 'pvalue' | 'genes'={} ]
    """

    cat2term2protein = {}
    for line in iterable:
         [Category, Term, Count, Percent,
         PValue, Genes, ListTotal,
         PopHits, PopTotal, FoldChange] = [x.strip for x in  line.split('\t')]
        if not d.has_key(Category):
            cat2term2protein[Category] = {}
        genes = Genes.strip(',').split(',')
        geneLs = {}
        for gene in genes: 
            geneLs[gene.strip().lstrip()] = True
        cat2term2protein[Category][Term] = {}
        cat2term2protein[Category][Term]['count'] = int(Count)
        cat2term2protein[Category][Term]['pvalue'] = float(PValue)
        cat2term2protein[Category][Term]['genes'] = geneLs
    return cat2term2protein

def cat2term2protein_sigFile(afile):
    """ Parse a DAVID chart file into
        category, term, proteins.
    @param afile: DAVID annotation output
    @return: [Category][Term]['count' | 'pvalue' | 'genes'={} ]
    """

    f_DAVID = open(afile)
    f_DAVID.readline()
    cat2term2protein = cat2term2protein_sig(f)
    f_DAVID.close()
    return cat2term2protein

def cat2term2protein_all(iterable):
    """ Take the output from DAVID annotation
        and organize by category, term, protein.
    
    @param iterable: open file or list that has DAVID output
    @return: [Category][Term]['count' | 'genes'={} ]
    """

    cat2term2protein = {}
    Categories = [x.strip for x in iterable[0].split('\t')[3:7]]
    for cat in Categories: cat2term2protein[cat] = {}
    for line in iterable[1:]:
        sp = [x.strip for x in line.split('\t')]
        for i in xrange(len(Categories)):
            cat = Categories[i]
            terms = [x.strip for x in sp[i+3].strip(',').split(':')]
            for x in xrange(len(terms)-1):
                if terms[x] != '':
                    term = terms[x].split()[-1]+':'+fullMonty(terms[x+1])
                    if not cat2term2protein[cat].has_key(term):
                        cat2term2protein[cat][term] = {}
                        cat2term2protein[cat][term]['count'] = 0
                        cat2term2protein[cat][term]['genes'] = {}
                    cat2term2protein[cat][term]['genes'][sp[0]] = True
                    cat2term2protein[cat][term]['count'] += 1               
    return cat2term2protein

def cat2term2protein_allFile(afile):
    """ Parse a DAVID annotation file into
        category, term, proteins.
    @param afile: DAVID annotation output
    @return: [Category][Term]['count' | 'genes'={} ]
    """

    f_DAVID = open(afile)
    cat2term2protein = cat2term2protein_all(f_DAVID)
    f_DAVID.close()
    return cat2term2protein

def getSigKEGGgenes_sigFile(afile, cutoff):
    """ Return {} of genes in significant KEGG
        pathways with p-values less than the
        cutoff, as described in the DAVID file.

    @param afile: DAVID chart file
    @param cutoff: p-value cut; less than #
    @return: {} of genes in sig KEGG pathways
    """

    genes = {}
    cat2term2protein = cat2term2protein_sigFile(afile)
    for path in cat2term2protein['KEGG_PATHWAY'].keys():
        if cat2term2protein['KEGG_PATHWAY'][path]['pvalue'] < cutoff:
            for g in cat2term2protein['KEGG_PATHWAY'][path]['genes'].keys():
                genes[g] = True
    return genes

def getSigKEGGgenes(motifSearchResults, bg, cutoff):
    """ Given a foreground and background {}, find
        significant KEGG pathways with DAVID and return
        {} of foreground genes that are in these pathways
        for the given p-value cutoff.
        
    @param motifSearchResults: foreground {}
    @param bg: background {}
    @param cuttoff: float p-value cutoff
    @return: {} of foreground genes in sig KEGG pathways
    """

    genes = {}
    cat2term2protein = getSigTerms(motifSearchResults, bg)
    for path in cat2term2protein['KEGG_PATHWAY'].keys():
        if cat2term2protein['KEGG_PATHWAY'][path]['pvalue'] < cutoff:
            for gene in cat2term2protein['KEGG_PATHWAY'][path]['genes'].keys():
                genes[gene] = True
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

def getAllTerms(geneFile):
    MozEmu = PyMozilla.MozillaEmulator(cacher=None,trycount=0)
    files = [['fileBrowser', 'testData', '']]
    fields = []
    MozEmu.download('http://david.abcc.ncifcrf.gov/tools.jsp')
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
    return cat2term2protein_all(MozEmu.download('http://david.abcc.ncifcrf.gov/UserDownload/' + addr  + '.txt').split('\n')[1:-1])

def getSigTerms(fgFile, bgFile):
    MozEmu = PyMozilla.MozillaEmulator(cacher=None,trycount=0)
    #try:
    files = [['fileBrowser', 'testData', '']]
    fields = []
    MozEmu.download('http://david.abcc.ncifcrf.gov/tools.jsp')
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

    userDownloadPage = MozEmu.download('http://david.abcc.ncifcrf.gov/chartReport.jsp?annot=,GOTERM_BP_5,GOTERM_CC_5,GOTERM_MF_5,KEGG_PATHWAY&currentList=0')
    index = userDownloadPage.find('UserDownload')
    addr = userDownloadPage[index+13:].split('.')[0]
    res = MozEmu.download('http://david.abcc.ncifcrf.gov/UserDownload/' + addr  + '.txt')
    return cat2term2protein_sig(res.split('\n')[1:-1])
