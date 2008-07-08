bimodal2pheno<-function(y,expData,threshold,bclambda,alpha){

#Find genes expressed in the "on" mode in a majority of samples 
#from a given phenotype

#Inputs
#y[n] = phenotype information 
#	y[i] = 1 if sample i is from the phenotype of interest, 0 otherwise 
#expData[d,n] = expression data, d = number of genes, n = number of samples
#threshold[d] = threshold separating low and high modes of expression
#bclambda[d] = box-cox lambda parameter
#alpha[d] = p-value for significance 

#Outputs
#onProbes = affymetrix probe set ids of genes expressed in the "on" mode
# in the majority of samples from a given phenotype 

#limit by phenotype 
expData<-expData[,y==1]

#convert threshold to expression space 
#reverse box-cox transformation 
threshold<-(threshold*bclambda+1)^(1/bclambda)

nclass<-sum(y)

#binomial test 
binaryData<-matrix(0,dim(expData)[1],dim(expData)[2])
binaryPval<-rep(NA,dim(expData)[1])
onIDX<-rep(0,dim(expData)[1])
for(i in 1:dim(expData)[1]){
	if(!is.nan(threshold[i])){
		binaryData[i,]<-expData[i,]>threshold[i]
		a=sum(binaryData[i,])
		binaryPval[i]<-1-pbinom(a,nclass,0.5)
		if(binaryPval[i]<=alpha) onIDX[i]<-1
	}
}

onProbes<-rownames(expData)[onIDX==1]

return(onProbes)
}
