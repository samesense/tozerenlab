nGap<-function(cl,data,clustMethod,B=500,k=2:25,method='complete'){

require(stats)

#Calculate the gap statistic to determine the optimal number of clusters 
#in a dataset

#Inputs
#cl = cluster object
#data[n,p] = expression data, n = number of samples, p = number of genes
#clustMethod: 1=kmeans
#	        2=hierarchical clustering  
#k = number of clusters for kmeans	

#Outputs 
#List with three components 
#gapStat[k] = gap statistic 
#sk[k] = standardized error of gap statistic
#khat = estimated optimum number of clusters

source('gap.stat.R')
clusterEvalQ(cl,library(stats))
clusterEvalQ(cl,source('gap.stat.R'))

clusterSetupRNG(cl,type='RNGstream')

partitions<-list()
if(clustMethod==1){
	clustFUN<-kmeans
	for(i in 1:length(k)){
		partitions[[i]]<-kmeans(data,k[i],nstart=10)$cluster
	}
}
if(clustMethod==2){
	clustFUN<-hclust
	tree<-hclust(dist(data),method=method)
	for(i in 1:length(k)){
		partitions[[i]]<-cutree(tree,(i+1))
	}
}

gap<-clusterApply(cl,partitions,gap.stat,data,B=B,clustFUN,clustMethod)

gapStat<-rep(NA,length(gap))
stError<-rep(NA,length(gap))
partitions<-matrix(NA,length(gap),nrow(data))
khat<-NULL
for(i in 1:length(gap)){
	gapStat[i]<-gap[[i]][1]
	stError[i]<-gap[[i]][2]
	partitions[i,]<-partitions[[i]]
	#identify khat
	if(i<length(gap) & is.null(khat)){
		if(gapStat[i]>=gap[[i+1]][1]-gap[[i+1]][2]) khat<-k[i]
	}
}

#partition<-partitions[[match(khat,k)]]
	
return(list(gap=gapStat,sk=stError,khat=khat))
}


	
