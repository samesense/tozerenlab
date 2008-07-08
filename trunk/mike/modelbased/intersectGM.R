intersectGM<-function(M,V,W){

#Find the intersection and overlap between component gaussians in 
#two-component gaussian mixture models

#Inputs
#M[d,2]=component means of Gaussian mixture models, d=number of genes
#V[d,2]=component standard deviations of Gaussian mixture models
#W[d,2]=component weights of Gaussian mixture models

#Outputs
#List with two components: 
#threshold[d] = point of intersection between component gaussians
#area[d,2] = area of intersection for each component gaussian

X<-rep(0,dim(M)[1])
A<-matrix(0,dim(M)[1],2)

#scale for accuracy
Scale<-M[,2]-M[,1]
M[,2]<-M[,1]+1
V<-V/matrix(Scale,dim(M)[1],2)

a<-V[,1]^2-V[,2]^2
b<-2*(M[,1]*V[,2]^2-M[,2]*V[,1]^2)
c<-(M[,2]*V[,1])^2-(M[,1]*V[,2])^2-2*(V[,1]*V[,2])^2*log((W[,2]*V[,1])/(W[,1]*V[,2]))

#roots
Xa<-matrix(0,dim(M)[1],2)
Xa[,1]<-((-b+(b^2-4*a*c)^(0.5))/(2*a))
Xa[,2]<-((-b-(b^2-4*a*c)^(0.5))/(2*a))

#check if root 1 >= mean of mode 1 and <= mean of mode 2
IDX<-Xa[,1]>=apply(M,1,min) & Xa[,1]<=apply(M,1,max)
IDX[is.na(IDX)]<-FALSE
X[IDX]<-Xa[IDX,1]

#check if root 2 >= mean of mode 1 and <= mean of mode 2
IDX<-Xa[,2]>=apply(M,1,min) & Xa[,2]<=apply(M,1,max)
IDX[is.na(IDX)]<-FALSE
X[IDX]<-Xa[IDX,2]

#area of overlap
A[,1] = 0.5*(erfc(abs(X - M[,1])/(V[,1]*2^0.5)))
A[,2] = 0.5*(erfc(abs(X - M[,2])/(V[,2]*2^0.5)))

#threshold
X = (X - M[,1])*Scale + M[,1]

return(list(threshold=X,area=A))
}

erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
