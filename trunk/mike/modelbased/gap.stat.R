gap.stat<-function(class,data,B,cluster.func,clustMethod){

#This function calculates the gap statistic to identify the number of clusters 
#present in the data.  It compares the partition provided by class to a reference
#null model of one cluster.
#Adapted from gap function in SAGx package

#data= nXp matrix of expression data, n=#of samples, p=# of genes 
#class= cluster indicator vector
#B= number of reference datasets 
#cluster.func= function used to cluster reference datasets 
#clustMethod=numeric indicator associated with cluster functions
#	1=kmeans
#	2=hierarchical clustering

require(stats)
FUN<-match.fun(cluster.func)

class.tab<-table(class)
nclus<-length(class.tab)
if(min(class.tab)==1) return(c(NA,NA))
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
	if(clustMethod==1){
		ref.clust<-cluster.func(z,nclus)
		ref.class<-ref.clust$cluster
	}
	if(clustMethod==2){
		ref.clust<-cluster.func(dist(z),method='complete')
		ref.class<-cutree(ref.clust,nclus)
	}
	
	tots[k]<-log(sum(by(z,factor(ref.class),intern<-function(x) sum(dist(x)/ncol(x))/2)))
}
out<-c(mean(tots)-temp1,sqrt(1+1/B)*sd(tots))
names(out)<-c('Gap statistic','one SE of simulation')
return(out)
} 
