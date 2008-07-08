wrap.multiDLDA<-function(cl,data,n.iterations,nfeatures,classes,Key=NULL,filename=NULL,random=F){

require(Biobase)
require(sma)
require(snow)
require(Rmpi)

#Supervised classification with multiple classes.  A separate binary 
#classifier is generated to identify each class versus all other samples.
#Test samples are classified according to the classifier with the highest 
#confidence 

#Inputs 
#cl = cluster object 
#data[n+1,p] = expression and phenotype data, n = number of samples, 
#	p = number of genes, first column equals phenotype information 	
#n.iterations = number of iterations of training and testing
#nfeatures = number of features to use for classification
#classes[nclasses] = 1:number of classes 
#Key[nclasses] = class descriptions, names(Key)=1:number of classes
#filename = filename used to output contingency table comparing prediction 
#	with true classification
#random = should features be selected randomly as control?

#Outputs 
#writes contingency table to filename 
#List with three components 
#max.class[nclasses,nclasses,length(nfeatures)]
#	 = contingency table comparing the true class to predicted class
#truth[nclasses,length(nfeatures)] = number of test samples in each class
#predictions[nclasses,length(nfeatures)] = number of predicted samples in 
#	each class 

clusterEvalQ(cl,library(Biobase))
clusterEvalQ(cl,library(sma))
clusterEvalQ(cl,source('multiDLDA.classification.R'))
clusterEvalQ(cl,source('multiRANDOM.R'))
clusterSetupRNG(cl, type='RNGstream')

n<-dim(data)[1]
p<-dim(data)[2]

#list of classes 
#each blade trains and tests on a different classifier
class.i<-list()
for(i in 1:length(classes)){
	class.i[[i]]<-classes[i]
}

#matrix of learning sets for each blade
#rows are samples, columns are iterations
y<-data[,1] 
ls.size<-round(2/3*n)
lset<-matrix(NA,nrow=ls.size,ncol=n.iterations)
for(i in 1:n.iterations){
	nclass<-nlevels(factor(y))
	props<-round(2/3*table(y))
	props[props==max(props)]<-max(props)-(sum(props)-ls.size)
	y.num<-as.numeric(factor(y))
	ls<-NULL
	for(j in 1:nclass){
		ls<-c(ls,sample(which(y.num==j))[1:props[j]])
	}
	lset[,i]<-ls
}
nclass.train<-props

nclass.test<-table(y[-ls])

if(random){
DLDA<-clusterApply(cl,class.i,multiRANDOM,data,nfeatures,lset)
}else{
DLDA<-clusterApply(cl,class.i,multiDLDA.classification,data,nfeatures,lset)
}

#tc=matrix(t-l,nit)
#pc=array(t-l,nit,nf)
#PIDs=matrix(nit,nf)

pc.matrix<-NULL
tc.matrix<-NULL
#determine predicted class
for(i in 1:length(DLDA)){ 
#loop over list elements 
	pc.i<-DLDA[[i]]$predicted.class 
	#convert to a vector 
	pc.matrix<-rbind(pc.matrix,as.vector(pc.i)) #rows,columns,height 

	tc.i<-DLDA[[i]]$y
	#convert to a vector 
	tc.matrix<-rbind(tc.matrix,as.vector(tc.i)) #rows,columns 
}

#test if rows of tc.matrix are equal
test.equal<-apply(tc.matrix,2,function(x) any(x!=x[1]))
if(any(test.equal)){
	stop('tc not equal between classifiers') 
}else{
	tc<-tc.matrix[1,]
}

#test for ties 
test.ties<-apply(pc.matrix,2,function(x) sum(x %in% max(x)))
if(any(test.ties>1)) print('ties')
#identify best classifier for each test sample, pc=vector length(rows*columns*height)
pc<-apply(pc.matrix,2,function(x) match(max(x),x))

#give max confidence level for each test sample 
pc.conf<-apply(pc.matrix,2,function(x) x[match(max(x),x)])

#divide vector into length(nfeatures) vectors
pc.nfeatures<-matrix(pc,nrow=length(nfeatures),ncol=length(pc)/length(nfeatures),byrow=T)

#nclass by nclass square matrix 
max.class<-array(0,dim=c(length(classes),length(classes),length(nfeatures)))
truth<-matrix(0,nrow=length(classes),ncol=length(nfeatures))
predictions<-matrix(0,nrow=length(classes),ncol=length(nfeatures))

for(i in 1:nrow(pc.nfeatures)){
	for(j in 1:ncol(pc.nfeatures)){
		temprow<-tc[j]
		tempcol<-pc.nfeatures[i,j]
		max.class[temprow,tempcol,i]<-max.class[temprow,tempcol,i]+1
		truth[temprow,i]<-truth[temprow,i]+1
		predictions[tempcol,i]<-predictions[tempcol,i]+1
	}
}

if(!is.null(Key)){
	for(i in 1:nrow(pc.nfeatures)){
		out<-max.class[,,i]/truth[,i]
		out<-round(out,2)
		out[out==0]<-'' #remove zeros from output 
		out<-rbind(Key,out)
		out<-cbind(c('',Key),out)
		filenames<-paste(filename,'_',nfeatures[i],'.txt',sep='')
		write(t(out),file=filenames,ncolumns=dim(out)[2],sep='\t')
	}
}
			
return(list(max.class=max.class,truth=truth,predictions=predictions))
}
