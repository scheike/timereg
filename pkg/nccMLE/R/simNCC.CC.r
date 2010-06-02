simNCC<-function(n,controls,lrr=c(0.25,-0.25),rate=1/0.004,win=0,d=1)
{
dimb<-length(lrr); 
zorig<-matrix(3*rnorm(dimb*n),n,dimb); 
time<-exp(zorig%*%lrr)*rexp(n)*rate;  
statuso<-(time<15); time[statuso==0]<-15

idx<-order(time); zorig<-as.matrix(zorig[idx,]); time<-time[idx]; statuso<-statuso[idx]; 
cases<-sum(statuso); 

samp<-c(); 
for (i in 1:cases)
{
is<-i+1;
id<-sample(is:n,size=controls,replace=F);
case<-c(i,time[i],zorig[i,],1); ctrl<-cbind(id,time[id],zorig[id,],0)
new<-rbind(case,ctrl); new<-cbind(new,i); samp<-rbind(samp,new); 
} 
#print(samp)

ud.f<-coxph(Surv(time,statuso) ~ zorig);
coxbeta<-ud.f$coef; varcox<-diag(ud.f$var); 
ud.ncc<-coxph(Surv(samp[,2],samp[,3+dimb]) ~ samp[,3:(3+dimb-1)]+strata(samp[,3+dimb+1]));
betancc<-ud.ncc$coef; varncc<-diag(ud.ncc$var); 
#print(ud.f);  print(ud.ncc)

caseid<-samp[samp[,3+dimb]==1,1]; contid<-samp[samp[,3+dimb]==0,1]; 
idtot<-unique(c(caseid,contid)); 
nobs<-length(idtot); 
Tobs<-c(time[1:cases],rep(9,nobs-cases));
status<-c(rep(1,cases),rep(0,nobs-cases)); 
nno<-n-nobs; Xobs<-zorig[idtot,]; 
m<-length(unique(contid)); 
cohort<-1*(!duplicated(c(unique(contid),idtot)));
cohort<-cohort[(m+1):(m+nobs)];

ud<-list(Xobs=Xobs,Tobs=Tobs,status=status,nno=nno)

return(ud); 
}
#ud<-simNCC(3728,2,rate=1/0.004,win=1)

simNCCg<-function(n,controls,lrr=-c(0.75,0.5,0.25,0.75,0.5,0.25),
rate=1/0.004,win=0,d=1,dimp=2,cor=0)
{
dimb<-length(lrr); zorig<-matrix(rnorm(dimp*n),n,dimp); 
vmat <- chol(matrix(c(1, cor, cor, 1), 2)); zorig<- zorig %*% vmat; 

zorig[,1]<-as.numeric(
cut(zorig[,1],breaks=quantile(zorig[,1]),include.lowest=T,labels=1:4)); 
zorig[,2]<-as.numeric(cut(zorig[,2],breaks=quantile(zorig[,2]),
include.lowest=T,labels=1:4)); 
zorig<-1*cbind(zorig[,1]==2, zorig[,1]==3, zorig[,1]==4,
zorig[,2]==2, zorig[,2]==3, zorig[,2]==4)
time<-exp(zorig%*%lrr)*rexp(n)*rate;  
statuso<-(time<15); time[statuso==0]<-15

idx<-order(time); zorig<-as.matrix(zorig[idx,]); time<-time[idx]; 
statuso<-statuso[idx]; cases<-sum(statuso); 

samp<-c(); 
for (i in 1:cases)
{
is<-i+1;
id<-sample(is:n,size=controls,replace=F);
case<-c(i,time[i],zorig[i,],1); ctrl<-cbind(id,time[id],zorig[id,],0)
new<-rbind(case,ctrl); new<-cbind(new,i); samp<-rbind(samp,new); 
} 
#print(samp)

ud.f<-coxph(Surv(time,statuso) ~ zorig);
coxbeta<-ud.f$coef; varcox<-diag(ud.f$var); 
ud.ncc<-coxph(Surv(samp[,2],samp[,3+dimb]) ~ samp[,3:(3+dimb-1)]+
strata(samp[,3+dimb+1]));
betancc<-ud.ncc$coef; varncc<-diag(ud.ncc$var); 
#print(ud.f);  print(ud.ncc)

caseid<-samp[samp[,3+dimb]==1,1]; contid<-samp[samp[,3+dimb]==0,1]; 
idtot<-unique(c(caseid,contid)); 
nobs<-length(idtot); 
Tobs<-c(time[1:cases],rep(9,nobs-cases));
status<-c(rep(1,cases),rep(0,nobs-cases)); 
nno<-n-nobs; Xobs<-zorig[idtot,]; 
m<-length(unique(contid)); 
cohort<-1*(!duplicated(c(unique(contid),idtot)));
cohort<-cohort[(m+1):(m+nobs)];

#print(Xobs); print(Tobs); print(status); print(nno); 

base0<-cbind(Tobs[status==1],(1/rate)*Tobs[status==1])
udem<-em.ncc(Xobs,Tobs,status,15,nno,Nit=100,beta=0,baseline=base0,
betait=5,emvar=1,d=d)
embeta<-udem$beta; emvar<-diag(udem$var); 

ud<-c(nobs,cases,coxbeta,varcox,betancc,varncc,embeta,emvar);
return(ud); 
}


simCC.ny<-function(n,m,lrr=-log(0.1),rate=1,r=0.5,tau=1,disk=1,win=0)
{
# n pop stoerrelse m cohort 
dimb<-length(lrr); 
if (disk==1) {zorig<-matrix((runif(dimb*n)>r)*1,n,dimb); 
for (j in 1:dimb) zorig[zorig[,j]==0,j]<-(-1);} else 
{zorig<-matrix(runif(dimb*n,-1,1),n,dimb);}; 
time<-exp(zorig%*%lrr)*rexp(n)*rate; statuso<-(time<tau); 
time[statuso==0]<-tau; 

idx<-order(time); zorig<-as.matrix(zorig[idx,]); time<-time[idx]; statuso<-statuso[idx]; 
cases<-sum(statuso); 

id.co<-sample(1:n,m,replace=F); id.case<-(1:n)[statuso==1];
z.case<-zorig[statuso==1,]; z.co<-zorig[id.co,];

id<-c(id.case,id.co); med<-!duplicated(id); id<-unique(id);
Xobs<-zorig[id,]; Tobs<-time[id]; status<-statuso[id]; Tno<-time[-id];
nno<-n-length(id);
cohort<-id; for (j in 1:length(id)) cohort[j]<-sum(id[j]==id.co);

base0<-cbind(Tobs[status==1],(1/rate)*Tobs[status==1])

ud<-list(Xobs=Xobs,Tobs=Tobs,status=status,nno=nno)

return(ud); 
} 


testII<-function(n,cc=F,ratei=1/0.004,stratsamp=F,controls=2,win=0) 
{
ratei<-200; n<-400; 
p<-2; 
stratumo<-sample(1:3,n,replace=T); 
zorig<-matrix(rnorm(p*n),n,p); 
xorig<-cbind(rbinom(n,1,0.5),sample(c(-10,0,10),n,replace=T)); 
obsrr<-exp(xorig%*% c(1,0.08))
time<-exp(zorig%*%c(0.24,-0.24))*rexp(n)*ratei*(1/obsrr); 
statuso<-(time<9); time[statuso==0]<-9; cases<-sum(statuso); 

idx<-order(time); zorig<-zorig[idx,]; time<-time[idx]; 
statuso<-statuso[idx]; stratumo<-stratumo[idx]; xorig<-xorig[idx]; 

print(table(stratumo,statuso))

m<-n*.2  # cohort size
if (cc==T) id.co<-sample(1:n,m,replace=F); 

if (cc==F) { ncc<-c(); 
if (stratsamp==F)  { 
for (i in 1:cases) {
is<-i+1; id<-sample(is:n,controls,replace=F); ncc<-c(ncc,i,id); } 

id.co<-!duplicated(c(1:cases,ncc)); mn<-length(id.co); 
id.co<-ncc[id.co[(cases+1):mn]]; } else { 
for (i in 1:cases) {
is<-i+1; 
id<-sample((is:n)[stratumo[is:n]==stratumo[i]],
controls,replace=F); ncc<-c(ncc,i,id); } 
id.co<-!duplicated(c(1:cases,ncc)); mn<-length(id.co); 
id.co<-ncc[id.co[(cases+1):mn]];  } 
}

id.case<-(1:cases); 
#z.case<-zorig[statuso==1,]; z.co<-zorig[id.co,]; 

id<-c(id.case,id.co); med<-!duplicated(id); id<-unique(id); 
Xobs<-as.matrix(zorig)[id,]; 
Tobs<-time[id]; status<-statuso[id]; Tno<-time[-id]; 
nobs<-length(id); nno<-n-nobs; 
stratum<-c(stratumo[id],stratumo[-id])

nrid<-1:n; 
#do<-cbind(1:n,stratumo,zorig,time,statuso)
#da<-cbind(c(nrid[id],nrid[-id]),stratum,c(zorig[id],zorig[-id]), 
#c(Tobs,Tno),c(status,statuso[-id]))
#cbind(da[order(da[,1]),],do)

ud<-list(Xobs=Xobs,Tobs=Tobs,status=status,nno=nno,stratum=stratum,
fullobsZ=xorig[-id])

XobsZ<-cbind(Xobs,xorig[id]); 

return(ud); 
}
