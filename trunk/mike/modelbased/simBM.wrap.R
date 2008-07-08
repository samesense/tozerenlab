simBM.wrap<-function(cl,n,p,iterations,nsg,nfeatures,mu,var1){

require(Biobase)
require(sma)
require(snow)
require(Rmpi)
require(igraph)

#create and classify simulated expression data

#Inputs
#cl = cluster object
#n = vector of number of observations/simulation
#p = number of features/simulation
#iterations = number of iterations/scenario
#nsg = vector of numbers of informative features
#nfeatures = vector of numbers of selected features
#mu = class-specific means
#var1 = class-specific variances

#Outputs 
#List of lists with four components 
#predicted.class[n,iterations] = score indicating predicted class,  
#predicted.genes[n,nfeatures,iterations] = genes identified as significant 
#y[n,iterations] = true classification
#true.genes[3,p,iterations] = first row equals true significant genes 
#	second row equals regression coefficient
#	third row equals variance 

clusterEvalQ(cl,library(Biobase))
clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(sma))
clusterEvalQ(cl,library(igraph))
clusterEvalQ(cl,source('simBM.classify.R'))
clusterSetupRNG(cl, type='RNGstream')

#Separate nfeatures and nsg for parallel processing
n.nsg<-matrix(NA,2,length(cl)) #matrix: [2,]=nsg,[1,]=n
n.nsg[2,]<-rep(nsg,(length(cl)/length(nsg))) #nsg[1],nsg[2],...
for(i in 1:ncol(n.nsg)){
	n.nsg[1,i]<-n[(i-1)%/%(length(cl)/length(n))+1] #rep(nsg[1],5),rep(nsg[2],5)
}
nsfsplit<-clusterSplit(cl,n.nsg)

nsim<-length(n)*length(nsg) #number of unique simulations
niterations<-iterations*nsim
cl.iterations<-niterations/length(cl)

classify<-clusterApply(cl,nsfsplit,simBM.classify,nfeatures,p,cl.iterations,mu,var1)
return(classify)
}

