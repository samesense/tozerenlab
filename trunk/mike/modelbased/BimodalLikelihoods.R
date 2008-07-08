BimodalLikelihoods<-function(expData,alpha){

require(mvtnorm)

#Identify bimodal genes using the likelihood ratio test. 
#Parameters of two component gaussian mixture models (GMM)
#estimated using expectation maximization

#Inputs
#expData[d,n] = box-cox transformed data,
# d = number of genes, n = number of arrays
#alpha = threshold indicating significance in likelihood ratio test 

#Output 
#list with six components: 
#W[d,2] = estimated weights of gaussian mixture model
#M[d,2] = estimated mean vectors of gaussian mixture model
#V[d,2] = estimated standard deviations of gaussian mixture model
#idx[d] = logical vector indicating bimodality

#parameters of components for Gaussian mixture model 
W<-matrix(0,nrow=dim(expData)[1],ncol=2)
M<-matrix(0,nrow=dim(expData)[1],ncol=2)
V<-matrix(0,nrow=dim(expData)[1],ncol=2)
L<-rep(0,dim(expData)[1])

#log-likelihoods
logL<-rep(0,dim(expData)[1])
logLGM<-rep(0,dim(expData)[1])

#p-value
pval<-rep(0,dim(expData)[1])
bimodalIDX<-rep(0,dim(expData)[1])

for(i in 1:dim(expData)[1]){
	#expectation-maximization
	bmParameters<-EM_GM(as.matrix(expData[i,]),k=2,ltol=0.05)
	tempM<-sort(bmParameters$M,index.return=T)
	W[i,tempM$ix]<-bmParameters$W
	M[i,]<-tempM$x
	V[i,tempM$ix[1]]<-sqrt(bmParameters$V[,,1])
	V[i,tempM$ix[2]]<-sqrt(bmParameters$V[,,2])
	L[i]<-bmParameters$L
		
	#sample mean and standard deviation
	mu<-mean(expData[i,])
	s<-sd(expData[i,])
	
	#log-likelihood 
	#unimodal
	normpdf<-dnorm(expData[i,],mean=mu,sd=s)
	logL[i]<-sum(log(normpdf[normpdf!=0]))

	#bimodal Gaussian mixture model
	pdfLM<-W[i,1]*dnorm(expData[i,],mean=M[i,1],sd=V[i,1])
	pdfHM<-W[i,2]*dnorm(expData[i,],mean=M[i,2],sd=V[i,2])
	pdfGM<-pdfLM+pdfHM
	logLGM[i]<-sum(log(pdfGM[pdfGM!=0]))

	#log-likelihood ratio test 
	testStat<--2*(logL[i]-logLGM[i])
	pval[i]<-1-pchisq(testStat,6)
	if(pval[i]<=alpha) bimodalIDX[i]<-1
}

return(list(W=W,M=M,V=V,idx=bimodalIDX))
}
 	
	