multiDLDA.classification<-function(class,data,nfeatures,lset){

#create classifier and classify test samples

#Inputs 
#class = number indicating class to build classifier for
#	classifiers for all classes are build using different parallel runs
#data[n+1,p] = expression and phneotype data
#nfeatures = number of features to use for classification
#lset[n*2/3,n.iterations] = learning sets to use for training 

#Outputs
#List of lists with four components 
#predicted.class[n,iterations] = score indicating predicted class,  
#predicted.genes[n,nfeatures,iterations] = genes identified as significant 
#y[n,iterations] = true classification
#phenotype = number indicating class to build classifier for

y2=data[,1] #class info
x1=data[,-(1)] #expression data
niterations<-dim(lset)[2]

#one vs. all 
y<-rep(0,length(y2))
y[y2==class]<-1

t.size<-length(y) #number of arrays
l.size<-round((2/3)*t.size) #size of learning set

true.class<-matrix(0,nrow=(t.size-l.size),ncol=niterations)
score<-array(0,dim=c((t.size-l.size),niterations,length(nfeatures))) #predicted classification
PIDs<-matrix(0,nrow=niterations,ncol=max(nfeatures)) #selected features

for(i in 1:niterations){
	ls<-lset[,i]
	true.class[,i]<-y2[-(ls)]
	
	#Select features by BWCSS/WCSS
	train.GE<-t(x1[ls,])
	rv1<-as.numeric(y[ls])
	rank.genes<-stat.bwss(train.GE,rv1) #rank genes by BWSS/WCSS
	rank.genes<-sort(rank.genes$bw,decreasing=T, index.return=T)
	PIDs[i,]<-rank.genes$ix[1:max(nfeatures)] #indices of top ranked genes

	for(m in 1:length(nfeatures)){
		#Classify data using diagonal linear discriminant analysis
		x<-t(x1[,PIDs[i,1:nfeatures[m]]])
		
		GELS<-x[,ls]
		GETS<-x[,-ls]
		rvLS<-as.matrix(y[ls])
	
		nls<-dim(GELS)[2]
		nts<-dim(GETS)[2]
		PDLS<-as.factor(rvLS)			
		index<-PDLS==levels(PDLS)[1] #0
		GE1<-GELS[,index]
		nls1<-dim(GE1)[2]
		index<-PDLS==levels(PDLS)[2] #1
		GE2<-GELS[,index]
		nls2<-dim(GE2)[2]

		#Calculate covariance
		GE1cov<-cov(t(GE1))
		GE2cov<-cov(t(GE2))
		GEcov<-((nls1-1)*GE1cov/(nls-2))+((nls2-1)*GE2cov/(nls-2))
		GEcov<-diag(GEcov)
		
		#Calculate distance to class centroids per gene
		dist<-matrix(0,nrow=nts, ncol=2)
		for(k in 1:nts){
			for(j in 1:length(GEcov)){
				GE1mean<-rowMeans(GE1)
				dist[k,1]<-dist[k,1]+((GETS[j,k]-GE1mean[j])^2)/GEcov[j]
				GE2mean<-rowMeans(GE2)
				dist[k,2]<-dist[k,2]+((GETS[j,k]-GE2mean[j])^2)/GEcov[j]
			}
		}
		
		#Make into distance formula
		score[,i,m]<-dist[,1]/(dist[,1]+dist[,2])
	}	
}
return(list(predicted.class=score,y=true.class,genes=PIDs,phenotype=class))
}
