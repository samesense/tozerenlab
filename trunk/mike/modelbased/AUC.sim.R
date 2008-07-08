AUC.sim<-function(sim.data){

#Calculate Area under the curve for simulated expression data 

#Inputs 
#sim.data = List of lists with four components 
#	predicted.class[n,iterations] = score indicating predicted class,  
#	predicted.genes[n,nfeatures,iterations] = genes identified as significant 
#	y[n,iterations] = true classification
#	true.genes[3,p,iterations] = first row equals true significant genes 
#	second row equals regression coefficient
#	third row equals variance 

#Output
#auc[length(sim.data)] = area under the curve for each simulation

y<-sim.data[[1]]$y
predicted.y<-sim.data[[1]]$predicted.class

require(ROC)
auc<-matrix(NA,ncol(y),length(sim.data))
step=0.01
spec.approx<-matrix(NA,(1/step+1),ncol(y))
sens.approx<-matrix(NA,(1/step+1),ncol(y))

for(j in 1:length(sim.data)){
	y<-sim.data[[j]]$y
	predicted.y<-sim.data[[j]]$predicted.class

	for(i in 1:ncol(y)){
	
	#Calculate specificity, sensitivity and AUC
	roc<-rocdemo.sca(y[,i],predicted.y[,i])
	auc[i,j]<-AUC(roc)

	spec<-1-roc@spec
	sens<-roc@sens
		if(spec[1]>spec[length(spec)]){
			spec<-rev(spec)
			sens<-rev(sens)
		}
	roc.approx<-approx(x=spec, y=sens, xout=seq(0,1,step))
	spec.approx[,i]<-roc.approx$x
	sens.approx[,i]<-roc.approx$y
	}
}
return(auc)
}
