GSA<-function(genelist,background,alpha){

#Uses hypergeometric test to assess enrichment of functional categories 
#Inputs
#genelist[d] = affymetrix probe ids for gene list of interest
#background[d,2] = first column equals affymetrix probe ids, second 
#	column equals probe annotation
#alpha = threshold indicating significance 

#Outputs
#sigterms[nterms,4] = first column equals significant terms, second
#	column equals observed number of probe sets in the term, 
#	third column equals expected number of probe sets in the term,
#	fourth column equals p-value derived from the hypergeometric test
#nterms = number of significant terms

#N = total number of affymetrix probes
#n = number of affymetrix probes in genelist 
#k = number of genes associated with current term 
#K = number of times each term appears 

N=dim(background)[1]

#find unique terms and number of times they appear in the reference list 
allterms<-NULL
for(t in 1:N){
	termtmp<-background[t,2]
	if(nchar(termtmp)!=0){
		multiterm<-unlist(gregexpr("///",termtmp))
		delimit<-c(-2, if(any(multiterm!=-1)) multiterm,nchar(termtmp)+1)
		for(D in 1:(length(delimit)-1)){
			allterms<-rbind(allterms,substr(termtmp,delimit[D]+3,delimit[D+1]-1))
		}
	}
}

allterms2<-unique(allterms)
print("Number of terms in reference list")
print(dim(allterms2)[1])
Klist<-rep(0,length(allterms2))
for(t in 1:dim(allterms2)[1]){
	Klist[t]<-sum(allterms==allterms2[t])
}

n<-length(genelist)
IDX<-match(genelist,background[,1])
listterms<-NULL

for(t in 1:length(IDX)){
	termtmp<-background[IDX[t],2]
	if(nchar(termtmp)!=0){
		multiterm<-unlist(gregexpr("///",termtmp))
		delimit<-c(-2,if(any(multiterm!=-1)) multiterm,nchar(termtmp)+1)
		for(D in 1:(length(delimit)-1)){
			listterms<-rbind(listterms,substr(termtmp,delimit[D]+3,delimit[D+1]-1))	
			#listgenes<-rbind(listgenes,c(" ",background[IDX[t],1]))
		}
	}
}

if(dim(listterms)[1]==0){
	print("No matches found")
	return()
}

listterms2<-unique(listterms)
print("Number of terms represented in genelist of interest")
print(dim(listterms2)[1])
K<-rep(0,dim(listterms2)[1])
E<-rep(0,dim(listterms2)[1])
O<-rep(0,dim(listterms2)[1])
P<-rep(0,dim(listterms2)[1])
#find number of times each term appears in the genelist 
for(t in 1:dim(listterms2)[1]){
	IDX<-listterms %in% listterms2[t]
	k<-sum(IDX)
	K[t]<-Klist[match(listterms2[t],allterms2)]
	E[t]<-n*K[t]/N
	O[t]<-k
	P[t]<-1-phyper((k-1),K[t],(N-K[t]),n)
}

#table of results 
sigterms<-NULL
for(t in 1:length(P)){
	if(P[t]<=alpha){
		sigterms<-rbind(sigterms,c(listterms2[t],O[t],E[t],P[t]))
	}	
}

return(sigterms)
}			
			