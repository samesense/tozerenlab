EM_GM<-function(X,k,ltol=0.1,maxiter=1000,pflag=0,Init=NULL){

#EM algortihm for k multidimensional Gaussian mixture estimation
#Adapted from code written for Matlab by 
#Patrick P. C. Tsui,
#PAMI research group
#Department of Electrical and Computer Engineering 
#University of Waterloo, 
#March, 2006 

#Inputs
#X[n,d] = input data, n=number of observations, d=dimension of variable
#k = maximum number of Gaussian components allowed
#ltol = percentage of the log likelihood difference between two iterations
#maxiter = maximum number of iteration allowed
#pflag = 1 for plotting GM for 1D or 2d cases only 
#Init = list of initial W,M,V: Init.W, Init.M, Init.V

#Outputs
#W[1,k] = estimated weights of GM
#M[d,k] = estimated mean vectors of GM
#V[d,d,k] = estimated covariance matrices of GM
#L = log-likelihood of estimates  

#Check input arguments
if(missing(X) || missing(k)){
	stop('EM_GM must have at least 2 inputs: X,k!')
}

#Initialize EM algorithm 
if(is.null(Init)){
	n<-dim(X)[1]
	d<-dim(X)[2]
	
	
	init.centers<-kmeans(X,centers=k,iter.max=100)
	Ci<-init.centers$cluster
	C<-init.centers$centers
	M<-t(C) #means

	Vp.count<-rep(0,k)
	Vp.X<-array(0,dim=c(n,d,k))
	for(i in 1:n){
		Vp.count[Ci[i]]<-Vp.count[Ci[i]]+1
		Vp.X[Vp.count[Ci[i]],,Ci[i]]<-X[i,]
	}
	V<-array(0,dim=c(d,d,k)) #variance
	W<-rep(0,k) #weights
	for(i in 1:k){
		W[i]<-Vp.count[i]/n;
		if(d==1){
			V[,,i]<-var(Vp.X[1:Vp.count[i],,i])
		}else{
			V[,,i]<-cov(Vp.X[1:Vp.count[i],,i])
		}	
	}
L<-0
}else{
	W<-Init$Init.W
	V<-Init$Init.V
	M<-Init$Init.M
}

#Initialize log likelihood 
Ln<-Likelihood(X,k,W,M,V)

Lo<-Ln*2

#### EM algorithm ####
niter<-0

while((abs(100*(Ln-Lo)/Lo)>ltol) & (niter<=maxiter)){
	#print(niter)
	E<-Expectation(X,k,W,M,V)
	Max<-Maximization(X,k,E)
	W<-Max$W
	M<-Max$M
	V<-Max$V
	Lo<-Ln
	Ln<-Likelihood(X,k,W,M,V)
	niter<-niter+1
}
L<-Ln
return(list(W=W,M=M,V=V,L=L))
}

Expectation<-function(X,k,W,M,V){
#This function is the modification of 'Expectation' in EM_GM made by Michael Boedigheimer to 
#enhance computational speed.  Note this modification takes more memory to execute. 
eps<-.Machine$double.eps
n<-dim(X)[1]
d<-dim(X)[2]
E<-matrix(0,n,k)
for(j in 1:k){
	if(all(V[,,j]==matrix(0,d,d))) V[,,j]<-matrix(1*eps,d,d)
	E[,j]<-W[j]*dmvnorm(X,M[,j],V[,,j])
}
total<-rep(rowSums(E),j)
E<-E/total
return(E)
}

Maximization<-function(X,k,E){
#This function is the modification of 'Maximization' in EM_GM made by Michael Boedigheimer to enhance
#computational speed. 
n<-dim(X)[1]
d<-dim(X)[2]
W<-colSums(E)
M<-t(X)%*%E/rep(W,d)
V<-array(0,dim=c(d,d,k))
for(i in 1:k){
	dXM<-X-rep(t(M[,i]),n)
	Wsp<-matrix(0,n,n)
	diag(Wsp)<-E[,i]
	V[,,i]<-t(dXM)%*%Wsp%*%dXM/W[i]
}
W<-W/n
return(list(W=W,M=M,V=V))
}


Likelihood<-function(X,k,W,M,V){
#Compute L based on K.V. Mardia, "Multivariate analysis", Academic Press, 1979, pg.96-97.
n<-dim(X)[1]
d<-dim(X)[2]
U<-t(mean(X))
S<-cov(X) #############################might need to fix######################

L<-0
for(i in 1:k){
	iV<-1/V[,,i] #################only works for d=1########################
	
	L=L+W[i]*(-0.5*n*log(2*pi*V[,,i])-0.5*(n-1)*(sum(diag(iV*S))+t(U-M[,i])*iV*(U-M[,i])))
}

return(L)
}
	




