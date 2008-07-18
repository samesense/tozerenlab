simBMNP<-function(n,p,iterations,nsg,nfeatures,mu,var1){

require(Biobase)
require(sma)
require(igraph)

#create and classify simulated expression data

#Inputs
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

classify.list<-list()
for(i in 1:length(n)){
	for(j in 1:length(nsg)){
		info<-c(n[i],nsg[j])
		classify<-simBM.classify(info,nfeatures,p,iterations,mu,var1)
		classify.list<-c(classify.list,list(classify))
	}
}

return(classify.list)
}

