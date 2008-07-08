simBM.classify<-function(nsf,nfeatures,p,iterations,mu,var1){

#create and classify simulated expression data 

#Inputs
#nsf[2] = c(number of simulated samples, number of significant genes)
#nfeatures = number of features to use for classification 
#p = number of features
#iterations = number of iterations
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

n<-nsf[1]
sg<-nsf[2]

sim.siggenes<-array(NA,dim=c(3,p,iterations))
true.class<-matrix(NA,nrow=n,ncol=iterations)
score<-matrix(NA,nrow=n, ncol=iterations)
PIDs<-array(0, dim=c(n,nfeatures, iterations))

mu1<-mu[1]
mu2<-mu[2]
j<-1
while(j<=iterations){
	
	#sim.data<-sim1(n,p,sg)
	#covariance structure: all genes independent
	covar<-matrix(0, nrow=p, ncol=p)
	diag(covar)<-rep(var1,p)
		
	#Select significant genes
	sig.genes<-sample(1:p, sg)
	y<-rbinom(n,1,0.5)
	
	#simulated expression data
	sim.data<-matrix(NA,n,p)
	sim.data[,-sig.genes]<-(1/2*(mvrnorm(p,rep(mu1,p),covar,empirical=T))+1/2*(mvrnorm(p,rep(mu2,p),covar,empirical=T)))[1:n,-sig.genes] #n*p matrix
	sim.data[y==0,sig.genes]<-mvrnorm(p,rep(mu1,p),covar,empirical=T)[sum(y==0),sig.genes]
	sim.data[y==1,sig.genes]<-mvrnorm(p,rep(mu2,p),covar,empirical=T)[sum(y==1),sig.genes]
 
	rownames(sim.data)<-1:n
	colnames(sim.data)<-1:p

	#Define regression coefficients, higher beta value associated with class 1
	beta<-rep(0,p)
	beta1<-runif(sg,0.1,1)
	sign<-sample(c(-1,1), sg, replace=T)
	beta1<-beta1*sign
	beta[sig.genes]<-beta1

	#Cross-multiply data with beta 
	terms<-as.vector(sim.data[,sig.genes] %*% beta[sig.genes]) #b1*x[,1] + b2*x[,2] + ...

	#logit(pi)=b0 + terms => pi=1/(1+e^(-b0 +terms))
	pi<-glm(y~1+offset(terms), family=binomial)$fitted.values
	
	sim.response<-rep(0,n)
	sim.response[pi>0.5]<-1
	
	if((sum(sim.response)>=(n-2))||(sum(sim.response)<=2)) next 

	genes<-matrix(0,nrow=3,ncol=p)
	genes[1,sig.genes]<-1
	genes[2,]<-beta
	genes[3,]<-diag(covar)
	colnames(genes)<-1:p
	rownames(genes)<-c('sig.genes','coef','covariance')
	
	#Output of sim
	x1<-sim.data
	y<-sim.response
	
	for(i in 1:nrow(x1)){
		ls<-(1:nrow(x1))[-i] #learning set	
		true.class[i,j]<-y[i] #class of test sample
		
		#feature selection by BWSS/WCSS
		train.GE<-t(x1[ls,])
		rv1<-as.numeric(y[-i])
		rank.genes<-stat.bwss(train.GE,rv1) #rank genes by BWSS/WCSS
		rank.genes<-sort(rank.genes$bw,decreasing=T, index.return=T)
	
		PIDs[i,,j]<-rank.genes$ix[1:nfeatures] #indices of top ranked genes
		x<-t(x1[,PIDs[i,,j]])	
		
		#Diagonal Linear Discriminant Analysis 
		GELS<-x[,ls]
		GETS<-x[,-ls]
		rvLS<-as.matrix(y[ls])

		nls<-dim(GELS)[2]
		PDLS<-as.factor(rvLS)
				
		index<-PDLS==levels(PDLS)[1]
		GE1<-GELS[,index]
		nls1<-dim(GE1)[2]
		
		index<-PDLS==levels(PDLS)[2]
		GE2<-GELS[,index]
		nls2<-dim(GE2)[2]

		#Calculate covariance
		GE1cov<-cov(t(GE1))
		GE2cov<-cov(t(GE2))
		GEcov<-((nls1-1)*GE1cov/(nls-2))+((nls2-1)*GE2cov/(nls-2))
		GEcov<-diag(GEcov)

		#Calculate distance to class centroids per gene
		dist<-rep(0,2)
		for(k in 1:length(GEcov)){
			GE1mean<-rowMeans(GE1)
			dist[1]<-dist[1]+((GETS[k]-GE1mean[k])^2)/GEcov[k]
			GE2mean<-rowMeans(GE2)
			dist[2]<-dist[2]+((GETS[k]-GE2mean[k])^2)/GEcov[k]
		}

		#Make into distance formula
  		score[i,j]<-dist[1]/(dist[1]+dist[2])
	}
	sim.siggenes[,,j]<-genes
	j<-j+1
}

return(list(predicted.class=score, predicted.genes=PIDs,  y=true.class, true.genes=sim.siggenes))
}



