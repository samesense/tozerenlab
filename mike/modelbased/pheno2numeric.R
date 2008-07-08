pheno2numeric<-function(phenotype){

#Convert phenotype information from string format to numeric

#Inputs
#phenotype[n,1] = phenotype information for all samples

#Outputs
#List with two components: 
#pheno[n] = phenotype in numeric format
#types[nphenotypes] = key linking number to string 
#nphenotypes = number of unique phenotypes 

types<-sort(unique(phenotype))

pheno<-rep(0,length(phenotype))
for(i in 1:length(types)){
	pheno[phenotype==types[i]]<-i
	names(types)[i]<-i
}

return(list(pheno=pheno,types=types))
}