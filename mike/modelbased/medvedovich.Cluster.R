medvedovich.Cluster<-function(out,groups,trueCi,nchains=25,max.iterations=NULL){

#Given the pairwise posterior distribution from all nodes find the average 
#pairwise posterior distribution, identify any outliers and cluster 

#Inputs
#out[niterations,8] = values of scan number, log-likelihood, K, alpha, lambda,
#	eta, a, b from every iteration of model based clustering	
#groups[niterations,n] = posterior samples of the group membership function
#trueCI[n] = phenotype information
#nchains = number of parallel chains used in execution 
#max.iterations = used to discard a number of iterations

#Outputs
#List with five components 
#PD[n,n] = pairwise posterior distribution 
#ci[n] = cluster indicator
#outliers = outlier samples
#ct = contingency table comparing cluster indicator to true classification
#ARI = adjusted rand index comparing cluster indicator to true classification

PD<-pairwisePD2(out,groups,nchains,max.iterations)$PD 
PD<-as.matrix(PD)

print('Percentage of posterior probabilities=0')
print(sum(PD==0)/(nsamples*nsamples))

#identify and remove outliers
outlierIDX<-rep(0,nsamples)
for(i in 1:nrow(PD)){
	if(sum(PD[i,]<=0.5)==(nsamples-1)){
		outlierIDX[i]<-1
	}
}
outliers<-(1:nsamples)[outlierIDX==1]
print('Outliers are')
print(outliers)
PD<-PD[!outlierIDX,!outlierIDX]
trueCi<-trueCi[!outlierIDX]

#distance measure 
D<-1-PD
D<-as.dist(D)

#hclust distances 
tree<-hclust(D,method='complete')
ci<-cutree(tree,h=0.99)
plot(tree)

#contingency table
cT<-table(trueCi,ci)
ARI<-adjustedRandIndex(trueCi,ci)
return(list(PD=PD,ci=ci,outliers=outliers,ct=cT,ARI=ARI))
}

pairwisePD2<-function(out,groups,nchains,max.iterations=NULL){

require(SparseM)

#Given a set of gibbs samplers run in parallel, find the most likely 
#partition from each chain and average them together

#Inputs
#out[niterations,8] = values of scan number, log-likelihood, K, alpha, lambda,
#	eta, a, b from every iteration of model based clustering	
#groups[niterations,n] = posterior samples of the group membership function
#trueCI[n] = phenotype information
#nchains = number of parallel chains used in execution 
#max.iterations = used to discard a number of iterations

#Outputs 
#List with three components 
#PD[n,n] = pairwise posterior distribution
#chainPD[nchains,n] = maximum likelihood partition from each run 
#chainLL = maximum likelihood from each run 

niterations<-dim(groups)[1]
nsamples<-dim(groups)[2]
niterChain<-niterations/nchains

if(is.null(max.iterations)){
	max.iterations<-niterChain
}else{
	max.iterations<-max.iterations
}

logL<-out[,2]
#find cluster from each partition with the max likelihood 
chainIDX<-1:max.iterations
maxLL<-rep(0,length(nchains))
clustChain<-matrix(0,nrow=nchains,ncol=nsamples)

for(i in 1:nchains){
	tempLL<-logL[chainIDX]
	tempGroups<-groups[chainIDX,]
	maxLL[i]<-tempLL[tempLL %in% max(tempLL)]
	clustChain[i,]<-tempGroups[tempLL %in% max(tempLL),]
	chainIDX<-chainIDX+niterChain
}

#create pairwise posterior distribution from max likelihood 
#partitions from each chain
tempPwise<-as.matrix.csr(0,nsamples,nsamples)
for(i in 1:nrow(clustChain)){
	pwiseChain<-as.matrix.csr(outer(clustChain[i,],clustChain[i,],'=='))
	tempPwise<-tempPwise+pwiseChain
}
pwisePD<-tempPwise/nrow(clustChain)

#outliers
if(any(apply(x,2,function(x) all(x<=0.5)))) print(outliers) 
	
return(list(PD=pwisePD,chainPD=clustChain,chainLL=maxLL))
}



