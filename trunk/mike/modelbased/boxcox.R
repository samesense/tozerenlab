boxcox<-function(data){

require(R.matlab)

#box-cox transformation, optimal lambda chosen using maximum likelihood 
#estimation run with fminsearch in matlab

#Inputs 
#data[d,n] = expression data matrix, d = number of genes, n = number of arrays

#Outputs
#list with two components:
#bct[d,n] = box-cox transformed expression data
#bclambda[1,d] = lambda parameter for box-cox transformation

#start matlab server
Matlab$startServer()
matlab<-Matlab()
Sys.sleep(10) #wait for server 
isOpen<-open(matlab)
on.exit(close(matlab))

#set data to Matlab variable 
setVariable(matlab,data=data)

#set function to boxcox transform data gene by gene 
setFunction(matlab,"                    				\
function [XFORM,LAMBDA]=bcXform(X)      				\
for i=1:size(X,1)				    				\
	[XFORM(i,:),LAMBDA(i)]=boxcox(double(X(i,:))');		\
end")

#get transformed data and lambda parameter 
evaluate(matlab,"[bct,bclambda]=bcXform(data);")
bcXform<-getVariable(matlab,c("bct","bclambda"))
return(bcXform)
}

