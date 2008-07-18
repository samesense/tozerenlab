nGapNP<-function(expData,clustMethod,B=500,k=2:25,method='complete'){

#Calculate the gap statistic to determine the optimal number of clusters 
#in a dataset

#Inputs
#expData[n,d] = expression data, n = number of samples, d = number of genes
#clustMethod = indicates which clustering method to use
#          1 = kmeans 
#          2 = hierarchical clustering
#k = number of clusters for kmeans 
#method = linkage method for hierarchical clustering 

#Outputs
#List with three components 
#gapStat[k] = gap statistic 
#sk[k] = standardized error of gap statistic
#khat = estimated optimum number of clusters

gap<-matrix(0,nrow=length(k),ncol=2)
for(i in 1:length(k)){
	
	if(clustMethod==1){
		clustFUN<-kmeans
		partition<-kmeans(expData,k[i],nstart=10)$cluster
	}
	if(clustMethod==2){
		clustFUN<-hclust 
		tree<-hclust(dist(expData),method=method)
		partition<-cutree(tree,k[i])
	}
	gap[i,]<-gap.stat(partition,expData,B=B,clustFUN,clustMethod)
}
print(gap)
khat<-NULL
for(i in 1:(dim(gap)[1])){
	if(i<length(gap) & is.null(khat)){
		if(is.na(gap[i+1,1])) next
		if(gap[i,1]>=gap[i+1,1]-gap[i+1,2]) khat<-k[i]
	}
	print(i)
}
print(gap)
stError=gap[,2]
gapStat=gap[,1]
return(list(gap=gapStat,sk=stError,khat=khat))
}
if(0){
gap.stat<-function(data,class,B=500,cluster.func){

#This function calculates the gap statistic to identify the number of clusters 
#present in the data.  It compares the partition provided by class to a reference
#null model of one cluster.
#Adapted from gap function in SAGx package

#Inputs
#data[n,p] = expression data, n = number of samples, p = number of genes 
#class[n] = cluster indicator 
#B = number of reference datasets 
#cluster.func = function used to cluster reference datasets 

#Outputs
#out[2] = c(gap.statistic, standardized error of statistic)

require(stats)
FUN<-match.fun(cluster.func)

class.tab<-table(class)
nclus<-length(class.tab)
if(min(class.tab)==1) stop('Singleton clusters not allowed')
if(!(length(class)==nrow(data))) stop('Length of class vector differs from nrow of data')

data<-as.matrix(data)
data<-scale(data,center=T,scale=F) #center genes so means==0
#generate reference distribution 
temp1<-log(sum(by(data,factor(class),intern<-function(x) sum(dist(x)/ncol(x))/2))) #log(Wk)
veigen<-svd(data)$v #singular value decomposition
x1<-crossprod(t(data),veigen) 

z1<-matrix(NA,nrow=nrow(x1),ncol=ncol(x1))
tots<-rep(NA,B)
for(k in 1:B){
	min.x<-apply(x1,2,min) #column mins
	max.x<-apply(x1,2,max) #column maxs
	z1<-matrix(runif(nrow(x1)*ncol(x1),min=min.x,max=max.x),nrow(x1),ncol(x1),byrow=T)
	z<-crossprod(t(z1),t(veigen)) #back-transform 
	ref.clust<-cluster.func(z,nclus)
	ref.class<-ref.clust$cluster
	tots[k]<-log(sum(by(z,factor(ref.class),intern<-function(x) sum(dist(x)/ncol(x))/2)))
}
out<-c(mean(tots)-temp1,sqrt(1+1/B)*sd(tots))
names(out)<-c('Gap statistic','one SE of simulation')
return(out)
} 

}