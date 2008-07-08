mbsc<-function(
                                              
 groups=ceiling(runif(dim(Y)[1],0,round(sqrt(dim(Y)[1])))),                       # starting values for groups, randomly assign to sqrt of N groups, for parallelization
Y,						# rows=objects, cols=attributes
 ps2=cbind(rep(1/2,dim(Y)[2]),apply(Y,2,var)/2),# gamma prior for base var
 pmu=cbind(apply(Y,2,mean),apply(Y,2,var)),     # normal prior for base mean
 peta=c(2,2),                                   # inv-gamma prior for eta
 ptheta=c(1,1),                                 # beta prior for pr relevance
 palpha=c(1,1),                                 # beta prior for alpha/alpha+1
 pwm=c(2,2),                                    # gamma prior for mean of 1/w^2
 pwv=c(2,2),                                    # gamma prior for var of 1/w^2
 nscan=1000,                                    # number of mcmc scans
 verb=T,                                        # verbose output?
 odens=max(1,round(nscan/1000)),                # how often to output results
 seed=1,                                        # seed of random numb gen 
 ngb=dim(Y)[1],                                 # no. gibbs samples per scan 
 nsm=5,                                         # no. of s-m proposals per scan 
 plt=F,                                         #
 ofile="CLUSTERS",
 saveres=T){                             # output filename

#### set seed
#set.seed(seed) #set by RNGstream

#set groups for convergence test 
if(length(groups)==1 & groups==1){
	groups=rep(1,dim(Y)[1])
}
if(length(groups)==1 & groups==2){
	groups=1:dim(Y)[1]
}

#### data size
n<<-dim(Y)[1] 
m<<-dim(Y)[2] 

#### starting values
mu<-pmu[,1]
s2<-pmu[,2]
E<-as.matrix(Y)
E<-t(t(E)-mu)
eta<-1
t2<-s2*eta
theta<- rep(log(ptheta[1]/ptheta[2]),m)
alpha<- palpha[1]/palpha[2]
a<-2 ; b<-2

#### Markov chain Monte Carlo
if(saveres==T) {
OUT<-matrix(NA,nrow=floor(nscan/odens),ncol=8) 
colnames(OUT)<-c("scan","ll","K","alpha","lambda","eta","miw","viw")
GROUPS<-matrix(NA,nrow=floor(nscan/odens),ncol=n)
MU<-S2<-matrix(NA,nrow=floor(nscan/odens),ncol=m)
uRm<-uDm<-uWm<-NULL
if(ngb+nsm==0) { uRt<-uDt<-uWt<-matrix(0,max(groups),m) }
              }
for(ns in 1:nscan) {

#### update marginal likelihood within each cluster
lpyk.marg<-rep(0,max(groups))
nk<-rep(0,max(groups))
for(k in 1:max(groups)){
nk[k]<-sum(groups==k)
Ek<-matrix(E[groups==k,],nrow=nk[k],ncol=m)
lpyk.marg[k]<-lpe(Ek,s2,t2,a,b,theta)
                       }
K<-k
####


#### update cluster memberships

## gibbs sampling
for(i in sample(1:n,ngb)){  
	
tmp<-rgroup.gibbs(i,nk,groups,K,lpyk.marg,E,s2,t2,a,b,theta,alpha) 
groups<-tmp$groups ; K<-tmp$K ; nk<-tmp$nk ; lpyk.marg<-tmp$lpyk.marg
                         }

## split merge
for(i in seq(1,nsm,length=nsm)){ 
tmp<-rgroups.sm.r(groups,E,s2,t2,a,b,theta,K,lpyk.marg,alpha)
groups<-tmp$groups;K<-tmp$K;nk<-tmp$nk;lpyk.marg<-tmp$lpyk.marg
                                }
####



#### update alpha
alpha<-r.alpha(n,K,palpha[1],palpha[2])
####

#### update parameters
RDW<-rrdw.clust(groups,E,s2,t2,a,b,theta)
R<-RDW$R ; W<-RDW$W*(R)+(1-R) ; WD<-RDW$D*R*W
                                                                                
## update mu
EE<- Y - WD
V<-  1/(apply(1/W^2,2,sum)/s2 + 1/pmu[,2])
M<-  V*( apply(EE/W^2,2,sum)/s2 + pmu[,1]/pmu[,2])
mu<-rnorm(m,M,sqrt(V))
E<-t(t(Y)-mu)
                                                                                
## update s2
EE<-t( t(Y-WD) -mu )
ss.s2<- apply( EE^2/W^2,2,sum) + apply( (RDW$uD)^2/eta,2,sum)
df.s2<- n + K
s2<-1/rgamma(m,ps2[,1]+df.s2/2,ps2[,2]+ss.s2/2)
####
                                                                                
#### update hyperparameters
                                                                                
## update theta
pr<-rbeta(1, ptheta[1]+sum(RDW$uR==1), ptheta[2]+sum(RDW$uR==0) )
theta<-rep(log(pr/(1-pr)),m)
                                                                                
## update t2/eta
eta<- 1/rgamma( 1, peta[1]+K*m/2, peta[2]+ sum(  t(t(RDW$uD)/sqrt(s2))^2)/2 )
t2<-s2*eta
                                                                                
## MH for a,b
miw<-a/b
viw<-a/b^2
                                                                                
w2<-c(RDW$uW)^2
em<-mean(1/w2)
vm<-var(1/w2)
nw<-m*K
miw.p<-rnorm(1,em,sqrt(vm/nw))
viw.p<-1/rgamma(1,1+nw/2,1+nw*vm/2 )
b.p<-miw.p/viw.p
a.p<-miw.p*b.p

lhr<- sum(dgamma(1/w2, a.p,b.p,log=T)) -  sum(dgamma(1/w2, a,b,log=T)) +
   dnorm(miw,em,sqrt(vm/nw),log=T) - dnorm(miw.p,em,sqrt(vm/nw),log=T) +
   dgamma(1/viw,1+nw/2,1+nw*vm/2,log=T)-dgamma(1/viw.p,1+nw/2,1+nw*vm/2,log=T)+
   dgamma(miw.p,pwm[1],pwm[2],log=T)-dgamma(miw,pwm[1],pwm[2],log=T)  +
   dgamma(viw.p,pwv[1],pwv[2],log=T)-dgamma(viw,pwv[1],pwv[2],log=T)
if(log(runif(1))<lhr) { a<-a.p ; b<-b.p }
####

#### output
if(ns%%odens==0){
out<-c(ns,sum(lpyk.marg),K,alpha,theta[1],eta,miw,viw)
cat(round( out,3),sort(table(groups)), "\n") 
if(saveres==T) {
OUT[ns/odens,]<-out
GROUPS[ns/odens,]<-groups
MU[ns/odens,]<-mu
S2[ns/odens,]<-s2
if(ngb+nsm==0) {
      uRt<-uRt+RDW$uR
      uDt<-uDt+RDW$uD
      uWt<-uWt+sqrt(RDW$uW)
      uRm<-uRt/( ns/odens) ; uDm<-uDt/( ns/odens) ; uWm<-uWt/( ns/odens)
                }
dscans<-1:(ns/odens)
tmp<-list(OUT=OUT[dscans,],GROUPS=GROUPS[dscans,],MU=MU[dscans,],S2=S2[dscans,],
          uRm=uRm,uDm=uDm,uWm=uWm)
dput(tmp,ofile) 
                }
if(plt==T){ 
 par(mfrow=c(3,2))
 for(i in c(2,3,5,6,7,8)){plot(OUT[,i],ylab=colnames(OUT)[i],typ="l")} 
           }

####
                }
               }
save.image(paste(ofile,".image",sep=""))   
list(OUT=OUT,GROUPS=GROUPS,MU=MU,S2=S2,uRm=uRm,uDm=uDm,uWm=uWm)
   }







########################################## Helper functions

system("R CMD SHLIB mbsc.c")

dyn.load("mbsc.so")

theta.hat<- function(Ek,s2,t2,a,b,theta) {
tmp<- .C("lpr",
      as.double(c(t(as.matrix(Ek)))),
      as.double(s2),
      as.double(t2),
      as.double(a),
      as.double(b),
      as.double(theta),
      as.integer(dim(Ek)[1]),
      as.integer(dim(Ek)[2]),
      th=double(dim(Ek)[2]),
      sm=double(1))
list(theta.hat=tmp$th,sm=tmp$sm)
               }

lpe<-function(Ek,s2,t2,a,b,theta) { theta.hat(Ek,s2,t2,a,b,theta)$sm }





#### sample S, Del given groupings
rrdw.clust<-function(groups,E,s2,t2,a,b,theta){
#this function needs rsdel.E
Rnew<-Dnew<-Wnew<-E*0 
uR<-uD<-uW<-NULL 

for(g in 1:max(groups)) {
i<-(1:n)[groups==g]
if(length(i)==1) { Eg<-matrix(E[i,],nrow=1,ncol=m) }
if(length(i)>1)  { Eg<-E[i,] }
tmp<-rrdw.E(Eg,s2,t2,a,b,theta)
uD<-rbind(uD,tmp$d) ; uR<-rbind(uR,tmp$r) ; uW<-rbind(uW,tmp$w)
Rnew[i,]<-matrix(tmp$r,nrow=length(i),ncol=m,byrow=T)
Dnew[i,]<-matrix(tmp$d,nrow=length(i),ncol=m,byrow=T)
Wnew[i,]<-matrix(tmp$w,nrow=length(i),ncol=m,byrow=T)
                  }
list(R=Rnew,D=Dnew,W=Wnew,uR=uR,uD=uD,uW=uW) }

#### helper for rrdw.clust
rrdw.E<-function(Eg,s2,t2,a,b,theta) {
ng<-length(c(Eg))/m

#sample r
theta.r<- theta.hat(Eg,s2,t2,a,b,theta)$theta.hat + theta
r<-rbinom(m,1,1/(1+exp(-theta.r)))

#sample w   
e1<-apply(Eg,2,mean)
e2<-apply(Eg^2,2,mean)
ngr<-ng*r
w<-1/sqrt(rgamma(m,a+ngr/2, b+.5*ngr*(e2-t2*e1^2/(s2/ng+t2))/s2 ))

#sample d
Ew<- t(t(Eg)/w)
vr<-1/(ngr/s2 + 1/t2)
mn<-vr*(  ngr*apply(Ew,2,mean)/s2 )
d<-rnorm(m,mn,sqrt(vr))
list( r=r,d=d,w=w) }

#### sample alpha
r.alpha<-function(n,K,a=1,b=1) {
th<-seq(0.001,.999,length=5000)
tmp<-K*log(th/(1-th))+lgamma( th/(1-th)) -lgamma(th/(1-th) +n) +
         dbeta(th,a,b,log=T) 
tmp<-exp(tmp-max(tmp))
rtheta<-sample(th,1,prob=tmp)
(rtheta)/(1-rtheta)
                        }


####


rgroup.gibbs<-function(i,nk,groups,K,lpyk.marg,E,s2,t2,a,b,theta,alpha ) {
# sample groups[i] conditional on  other values
nk[groups[i]]<-nk[groups[i]]-1
al<-rep(0,K)

if(sum(groups==groups[i])==1) {
lpy.mi<-c(lpyk.marg)
lpy.mi[groups[i]]<-0
lpy.pi<-rep(NA,K)
lpy.pi[groups[i]]<-lpyk.marg[groups[i]]

for(k in (1:K)[-groups[i]]) {
      i.k<-(1:n)[groups==k]
      lpy.pi[k]<-lpe(rbind(E[i.k,],E[i,]),s2,t2,a,b,theta)
                             }
al[groups[i]]<-alpha           }

if(sum(groups==groups[i])>1) {
lpy.mi<-c(lpyk.marg)
lpy.pi<-rep(NA,K)
lpy.pi[groups[i]]<-lpyk.marg[groups[i]]

ik.mi<-(1:n)[groups==groups[i]]
ik.mi<-ik.mi[ik.mi!=i]
lpy.mi[groups[i]]<-lpe(E[ik.mi,,drop=F],s2,t2,a,b,theta)

for(k in (1:K)[-groups[i]]) {
      i.k<-(1:n)[groups==k]
 lpy.pi[k]<-lpe(rbind(E[i.k,],E[i,]),s2,t2,a,b,theta)
                             }
lpy.mi<-c(lpy.mi,0)
lpy.pi<-c(lpy.pi,lpe(E[i,,drop=F],s2,t2,a,b,theta) )
K<-K+1
nk<-c(nk,0)
al<-c(al,alpha)            }

pki<- lpy.pi-lpy.mi + log( nk+al )

pki<- exp( pki-max(pki) )/sum( exp( pki-max(pki) )  )
ki<-sample( 1:K, 1, prob=pki)

lpyk.marg<-lpy.mi*( (1:K)!=ki) + lpy.pi*( (1:K)==ki)
groups[i]<-ki
nk[ki]<-nk[ki]+1
lpyk.marg<-lpyk.marg[unique(groups)]
nk<-nk[unique(groups)]
groups<-match(groups,unique(groups))
K<-max(groups)
list(groups=groups,K=max(groups),nk=nk, lpyk.marg=lpyk.marg)
 }


#### split-merge proposals
rgroups.sm.r<-function(groups,E,s2,t2,a,b,theta,K,lpyk.marg,alpha) {


ij<-sample(1:n,2)
i<-ij[1]
j<-ij[2]
                                                                                
lhr<-0
                                                                                
if(groups[i]==groups[j]) { #split it
split<-T
groups.p<-groups
groups.p[i]<-max(groups)+1
groups.p[j]<-max(groups)+2
tmp<-NULL
lset<-((1:n)[groups==groups[i] ] )
lset<-lset[lset!=i & lset!=j]
groups.p[lset]<-0
if(length(lset>0)){
for( l  in sample(lset)) {
                                                                                
lci<-(1:n)[groups.p==groups.p[i]]
lcj<-(1:n)[groups.p==groups.p[j]]
lpr<- (  lpe(E[c(lcj,l),],s2,t2,a,b,theta) -
         lpe(E[c(lcj),,drop=F],s2,t2,a,b,theta)  ) -
      (  lpe(E[c(lci,l),],s2,t2,a,b,theta) -
         lpe(E[c(lci),,drop=F],s2,t2,a,b,theta)  )
nr<-sum(groups.p==groups.p[j])/sum(groups.p==groups.p[i])
pil<-1/(1+  nr* exp(lpr) )
lg<-rbinom( 1,1,pil)
groups.p[ l]<-groups.p[i]*lg + groups.p[j]*(1-lg)
lhr<-lhr + log(c(1-pil,pil))[lg+1]
                           }
                       }
                          }
if(groups[i]!=groups[j]) {
split<-F
groups.o<-groups.p<-groups
groups.p[groups==groups[j]]<-groups[i]
                                                                                
tmp<-NULL
lset<-((1:n)[groups.p==groups.p[i] ] )
lset<-lset[lset!=i & lset!=j]
groups.o[lset]<-0
if(length(lset>0)){
for( l  in sample(lset)) {
lci<-(1:n)[groups.o==groups.o[i]]
lcj<-(1:n)[groups.o==groups.o[j]]
                                                                                
lpr<- (  lpe(E[c(lcj,l),],s2,t2,a,b,theta) -
         lpe(E[c(lcj),,drop=F],s2,t2,a,b,theta)  ) -
      (  lpe(E[c(lci,l),],s2,t2,a,b,theta) -
         lpe(E[c(lci),,drop=F],s2,t2,a,b,theta)  )
nr<-sum(groups.o==groups.o[j])/sum(groups.o==groups.o[i])
                                                                                
pil<-1/(1+  nr* exp(lpr) )
lg<-1*( groups[l]==groups[i] )
groups.o[l]<-groups[l]
lhr<-lhr + log(c(1-pil,pil))[lg+1]
                           }
                      }
                          }
##########
groups.p<-match(groups.p,unique(groups.p))
K.p<-max(groups.p)
lpyk.marg.p<-nk<-rep(0,K.p)
for(k in 1:K.p){
nk[k]<-sum(groups.p==k)
lpyk.marg.p[k]<-lpe( E[groups.p==k,,drop=F],s2,t2,a,b,theta)
                       }

list(groups=groups,K=max(groups),nk=table(groups),lpyk.marg=lpyk.marg)
}                     




nmatch<-function(groups.0,groups.1) {
n<-length(groups.0)
G0<-matrix(groups.0,nrow=n,ncol=n,byrow=T)
M0<-1*(G0==t(G0)) ; diag(M0)<-NA
G1<-matrix(groups.1,nrow=n,ncol=n,byrow=T)
M1<-1*(G1==t(G1)) ; diag(M1)<-NA
 
 
Ndd<-sum( M0==0 & M1==0,na.rm=T)/2
Nsd<-sum( M0==1 & M1==0,na.rm=T)/2
Nds<-sum( M0==0 & M1==1,na.rm=T)/2
Nss<-sum( M0==1 & M1==1,na.rm=T)/2
list(Ndd=Ndd,Nsd=Nsd,Nds=Nds,Nss=Nss)                  }
 
 
jac.index<-function(groups.0,groups.1) {
tmp<-nmatch(groups.0,groups.1)
tmp$Nss/(tmp$Nds+tmp$Nsd+tmp$Nss)     }
 

add.output2<-function(){
par(mfrow=c(3,2))
for( g  in (1:(max(groups)))[table(groups)>1]   ) {
Yg<-Y[groups==g,,drop=F]
plot( c(1,m),range( t(Y)-mu),type="n")
mtext( dim(Yg)[1] ,side=3)
for(i in 1:dim(Yg)[1] ) {lines(1:m,Yg[i,]-mu,col="yellow") }
lines( 1:m, apply( t(t(Yg)-mu),2,mean),col="brown")
lines( 1:m, (S*Del)[min((1:n)[groups==g]),],col="blue")
                                    }
                          }

add.output<-function(){
par(mfrow=c(3,2))
plot(OUT[,2], type="l") 
plot(OUT[,3], type="l") 
plot(OUT[,5], type="l") 
plot(OUT[,6], type="l") 
plot(OUT[,7], type="l") 
plot(OUT[,8], type="l") 
                                                                                
                         }



 
nmatch<-function(groups.0,groups.1) {
n<-length(groups.0)
G0<-matrix(groups.0,nrow=n,ncol=n,byrow=T)
M0<-1*(G0==t(G0)) ; diag(M0)<-NA
G1<-matrix(groups.1,nrow=n,ncol=n,byrow=T)
M1<-1*(G1==t(G1)) ; diag(M1)<-NA
 
 
Ndd<-sum( M0==0 & M1==0,na.rm=T)/2
Nsd<-sum( M0==1 & M1==0,na.rm=T)/2
Nds<-sum( M0==0 & M1==1,na.rm=T)/2
Nss<-sum( M0==1 & M1==1,na.rm=T)/2
list(Ndd=Ndd,Nsd=Nsd,Nds=Nds,Nss=Nss)                  }
 
 
jac.index<-function(groups.0,groups.1) {
tmp<-nmatch(groups.0,groups.1)
tmp$Nss/(tmp$Nds+tmp$Nsd+tmp$Nss)     }
 

groups.mode<-function(GROUPS,nm=10) {
n<-dim(GROUPS)[2]
gcode<-runif(n)
gid<- GROUPS%*%gcode
tmp<-table(gid)
ttmodes<-(-sort(-tmp))[1:10]
GROUPS.modes<-NULL
for(i in 1:10) {
mid<-round(as.numeric(names(ttmodes)[i]),4)
GROUPS.modes<-rbind(GROUPS.modes,(GROUPS[ round(gid,4)==mid,,drop=F])[1,] )
                }
p.modes<-as.vector(ttmodes/dim(GROUPS)[1] )
list( GROUPS.modes=GROUPS.modes, p.modes=p.modes )
                                     }

groups.mean<-function(GROUPS) {

Dk<-matrix(0,nrow=dim(GROUPS)[1],ncol=dim(GROUPS)[1])
for(i in 1:(dim(Dk)[1]-1)) {
for(j in (i+1):dim(Dk)[1]) {
     Dk[i,j]<-Dk[j,i]<- Kdist(GROUPS[i,],GROUPS[j,])
                           }}
tdist<-apply(Dk^2,1,mean)
i.mean<-((1:dim(GROUPS)[1])[tdist==min(tdist)])[1]
groups.mean<-GROUPS[i.mean,]
list(Dist=Dk, groups.mean=groups.mean)
                                }

Kdist<-function(groups.0,groups.1) {
n<-length(groups.0)
G0<-matrix(groups.0,nrow=n,ncol=n,byrow=T)
M0<-1*(G0==t(G0)) ; diag(M0)<-NA
G1<-matrix(groups.1,nrow=n,ncol=n,byrow=T)
M1<-1*(G1==t(G1)) ; diag(M1)<-NA
                                                                                
                                                                                
Nsd<-sum( M0==1 & M1==0,na.rm=T)/2
Nds<-sum( M0==0 & M1==1,na.rm=T)/2
#2*(Nsd+Nds)
(Nsd+Nds)/(choose(n,2))
                     }


groups.conc<-function(GROUPS,cutoff=.5) {
n<-dim(GROUPS)[2]
MG<-matrix(0,n,n)
for(l in 1:dim(GROUPS)[1]){
g<-GROUPS[l,]
G<-matrix(g,nrow=n,ncol=n,byrow=T)
M<-1*(G==t(G)) ; diag(M)<-NA
MG<-MG+M         }
MG<-MG/(dim(GROUPS)[1])
                                                                                
#now complete the graph
CG<- 1*(MG>cutoff)
                                                                                
groups.c<-1
for(i in 2:n) {
mgi<-MG[i,1:(i-1)]
gp<-(1:i)[mgi>.5]
if(length(gp)>0) { groups.c<-c(groups.c, groups.c[gp[1]]) }
if(length(gp)==0) { groups.c<-c(groups.c, max(groups.c)+1) }
                }
groups.c }

sort.cluster<-function(groups) {match( groups,order(-table(groups))) }

