aalenBaseTest<-function(times,fdata,designX,status,id,clusters,
                        robust=0,sim=0,retur=0,antsim=1000,
                        weighted.test=1,covariance=0,resample.iid=0,
                        namesX=NULL,offsets=0,weights=0,silent=0){
  Ntimes<-length(times)
  designX<-as.matrix(designX); 
  if(is.matrix(designX) == TRUE) p <- as.integer(dim(designX)[2])
  if(is.matrix(designX) == TRUE) nx <- as.integer(dim(designX)[1])

  if (robust==0 & sim>=1)  robust<-1;  devi<-rep(0,3); 
  cumint<-matrix(0,Ntimes,p+1); Vcumint<-cumint; 
  if (retur==1) cumAi<-matrix(0,Ntimes,fdata$antpers) else cumAi<-0;
  test<-matrix(0,antsim,3*p); testOBS<-rep(0,3*p); testval<-c(); 
  unifCI<-c();
  rani<--round(runif(1)*10000)


                                        # 50 tilfaeldige score processer til test H: b(t)=b returneres
  if (sim>=1) simUt<-matrix(0,Ntimes,50*p) else simUt<-NULL;
  Ut<-matrix(0,Ntimes,p+1); 
  if (covariance==1) covs<-matrix(0,Ntimes,p*p) else covs<-0;
  if (resample.iid==1) {
    B.iid<-matrix(0,Ntimes,fdata$antclust*p) } else B.iid<-NULL;
  if (sum(offsets)==0) mof<-0 else mof<-1;
  if (sum(weights)==0) mw<-0 else mw<-1;

  robVar<-Vcumint;
  aalenout<- .C("robaalentest",
                as.double(times), as.integer(Ntimes),
                as.double(designX),as.integer(nx),as.integer(p),
                as.integer(fdata$antpers),as.double(fdata$start),as.double(fdata$stop),
                as.double(cumint),as.double(Vcumint),as.double(robVar),
                as.integer(sim),as.integer(antsim),as.integer(retur),
                as.double(cumAi),as.double(test),as.integer(rani),as.double(testOBS),
                as.integer(status),as.double(Ut),as.double(simUt),as.integer(id),
                as.integer(weighted.test),as.integer(robust),as.integer(covariance),
                as.double(covs),as.integer(resample.iid),as.double(B.iid),
                as.integer(clusters),as.integer(fdata$antclust),as.double(devi),
                as.integer(mof),as.double(offsets),as.integer(mw),as.double(weights),
                as.integer(silent)
                ,PACKAGE="timereg")

  if (covariance==1)  {
    covit<-matrix(aalenout[[26]],Ntimes,p*p);
    cov.list<-list();
    for (i in 1:Ntimes) cov.list[[i]]<-matrix(covit[i,],p,p);
  } else cov.list<-NULL;

  if (resample.iid==1)  {
    covit<-matrix(aalenout[[28]],Ntimes,fdata$antclust*p); 
    B.iid<-list(); 
    for (i in (0:(fdata$antclust-1))*p) {
      B.iid[[i/p+1]]<-as.matrix(covit[,i+(1:p)]); 
      colnames(B.iid[[i/p+1]])<-namesX; } 
  }

  robV<-matrix(aalenout[[11]],Ntimes,p+1);
  if (retur==1) {
    cumAi<-matrix(aalenout[[15]],Ntimes,fdata$antpers*1); 
    cumAi<-list(time=times,dM=cumAi,dM.iid=cumAi); 
                                        #cumAI<-list(); 
                                        #for (i in (0:(fdata$antpers-1))*p) cumAI[[i/p+1]]<-cumAi[,i+(1:p)]
  } else cumAi<-NULL;

  if (sim>=1) {
    Uit<-matrix(aalenout[[21]],Ntimes,50*p); 
    UIt<-list(); for (i in (0:49)*p) 
      UIt[[i/p+1]]<-as.matrix(Uit[,i+(1:p)]);
    Ut<-matrix(aalenout[[20]],Ntimes,(p+1)); 
    test<-matrix(aalenout[[16]],antsim,3*p);
    testOBS<-aalenout[[18]];
    for (i in 1:(3*p)) testval<-c(testval,pval(test[,i],testOBS[i]))
    for (i in 1:p) unifCI<-as.vector(c(unifCI,percen(test[,i],0.95)));
    pval.testBeq0<-as.vector(testval[1:p]);
    pval.testBeqC<-as.vector(testval[(p+1):(2*p)]);
    pval.testBeqC.is<-as.vector(testval[(2*p+1):(3*p)]);
    obs.testBeq0<-as.vector(testOBS[1:p]);
    obs.testBeqC<-as.vector(testOBS[(p+1):(2*p)]);
    obs.testBeqC.is<-as.vector(testOBS[(2*p+1):(3*p)]);
    sim.testBeq0<-as.matrix(test[,1:p]);
    sim.testBeqC<-as.matrix(test[,(p+1):(2*p)]);
    sim.testBeqC.is<-as.matrix(test[,(2*p+1):(3*p)]);
  } else {test<-NULL; unifCI<-NULL; Ut<-NULL; UIt<-NULL; 
          pval.testBeq0<-NULL;pval.testBeqC<-NULL; obs.testBeq0<-NULL;obs.testBeqC<-NULL;
          sim.testBeq0<-NULL;sim.testBeqC<-NULL; 
          sim.testBeqC.is<-NULL; pval.testBeqC.is<-NULL; 
          obs.testBeqC.is<-NULL; 
        }
  cumint <-matrix(aalenout[[9]],Ntimes,p+1);
  Vcumint<-matrix(aalenout[[10]],Ntimes,p+1);
  devi<-aalenout[[31]]

  list(cum=cumint,var.cum=Vcumint,robvar.cum=robV,residuals=cumAi,
       pval.testBeq0=pval.testBeq0, obs.testBeq0=obs.testBeq0,
       pval.testBeqC=pval.testBeqC, pval.testBeqC.is=pval.testBeqC.is,
       obs.testBeqC=obs.testBeqC,obs.testBeqC.is=obs.testBeqC.is,
       sim.testBeq0= sim.testBeq0,
       sim.testBeqC=sim.testBeqC,sim.testBeqC.is=sim.testBeqC.is,
       conf.band=unifCI,test.procBeqC=Ut,sim.test.procBeqC=UIt,
       covariance=cov.list,B.iid=B.iid,deviance=devi)
}

semiaalenTest<-function(times,fdata,designX,designG,status,id,clusters,
bhat=0,gamma=0,robust=0,sim=0,antsim=1000,weighted.test=1,retur=0,
covariance=0,
resample.iid=0,namesX=NULL,namesZ=NULL,deltaweight=1,weights=0,offsets=0,
fixgam=0,pseudo.score=0,silent=0)
{
Nalltimes <- length(times);  
Ntimes<-sum(status[(fdata$stop>times[1]) & (fdata$stop<=times[Nalltimes])])+1;
#print(Ntimes); print(Nalltimes); 
#print(times);  print(status);  print(cbind(fdata$stop,fdata$start,status))
designX<-as.matrix(designX); designG<-as.matrix(designG); 
if(is.matrix(designX) == TRUE) px <- as.integer(dim(designX)[2])
if(is.matrix(designX) == TRUE) nx <- as.integer(dim(designX)[1])
if(is.matrix(designG) == TRUE) pg <- as.integer(dim(designG)[2])
if(is.matrix(designG) == TRUE) ng <- as.integer(dim(designG)[1])
nb<-1; 
if(is.matrix(bhat)==TRUE) nb<-as.integer(dim(bhat)[1]); 
if(is.matrix(bhat)==FALSE) bhat<-matrix(0,nb,px+1); 
if (sum(offsets)==0) mof<-0 else mof<-1;
if (sum(weights)==0) mw<-0 else mw<-1;

if (covariance==1) covs<-matrix(0,Ntimes,px*px) else covs<-0;

if (pseudo.score>=1) { 
pscoret<-matrix(0,Ntimes,pg+1);
pscoreiid<-matrix(0,Ntimes,fdata$antclust*pg);
IintZHZt<-matrix(0,Ntimes,pg*pg); ptest<-matrix(0,pseudo.score,pg);} 
else {pscoret<-pscoreiid<-IintZHZt<-ptest<-0}

if (resample.iid==1) {
gamma.iid<-matrix(0,fdata$antclust,pg); 
B.iid<-matrix(0,Ntimes,fdata$antclust*px) }
else { B.iid<-gamma.iid<-NULL; }

if (retur==1) cumAi<-matrix(0,Ntimes,fdata$antpers) else cumAi<-0;
if (nx!=ng) print(" A design og B designs er ikke ens\n");
cum<-matrix(0,Ntimes,px+1); 
Vcum<-cum; robVcum<-cum; 
if (sum(abs(gamma))==0) gamma<-rep(0,pg)  else gamma<-gamma; 
Vargam<-matrix(0,pg,pg); RobVargam<-matrix(0,pg,pg); 
intZHZ<-matrix(0,pg,pg); intZHdN<-rep(0,pg); 

test<-matrix(0,antsim,3*px); testOBS<-rep(0,3*px); testval<-c(); 
rani<--round(runif(1)*10000); devi<-rep(0,3); 

# 50 tilfaeldige score processer til test H: b(t)=b returneres
if (sim==1) simUt<-matrix(0,Ntimes,50*px) else simUt<-NULL;
Ut<-matrix(0,Ntimes,px+1); 
antpers<-fdata$antpers;
intHdN<-rep(0,antpers); intHZ<-matrix(0,antpers,pg); 
varintHdN<-matrix(0,antpers,antpers); 

semiout<-.C("semiaalentest",
as.double(times),as.integer(Nalltimes),as.integer(Ntimes),
as.double(designX),as.integer(nx),as.integer(px),
as.double(designG),as.integer(ng),as.integer(pg),
as.integer(antpers),as.double(fdata$start),as.double(fdata$stop), 
as.integer(nb),as.double(bhat),
as.double(cum),as.double(Vcum),as.double(robVcum),
as.double(gamma),as.double(Vargam),as.double(RobVargam),
as.integer(sim),as.integer(antsim), 
as.double(test),as.integer(rani),as.double(testOBS),as.integer(robust),
as.integer(status),as.double(Ut),as.double(simUt),as.integer(id),
as.integer(weighted.test),as.double(cumAi),as.integer(retur),as.integer(covariance),as.double(covs),
as.integer(resample.iid),as.double(gamma.iid),as.double(B.iid),
as.integer(clusters),as.integer(fdata$antclust),as.double(devi),
as.double(intZHZ),as.double(intZHdN),as.integer(deltaweight), 
as.integer(mof),as.double(offsets),
as.integer(mw),as.double(weights),
as.integer(fixgam),
as.integer(pseudo.score),as.double(pscoret),as.double(pscoreiid),
as.double(IintZHZt),as.double(ptest),
as.double(intHdN),as.double(intHZ),as.double(varintHdN),as.integer(silent),
PACKAGE="timereg"); 

if (resample.iid==1)  {
gamma.iid<-matrix(semiout[[37]],fdata$antclust,pg);
covit<-matrix(semiout[[38]],Ntimes,fdata$antclust*px); 
B.iid<-list(); 
for (i in (0:(fdata$antclust-1))*px) {
B.iid[[(i/px)+1]]<-as.matrix(covit[,i+(1:px)]);
colnames(B.iid[[i/px+1]])<-namesX; }
}
devi<-semiout[[41]]

if (pseudo.score>=1)  {
IintZHZt<-matrix(semiout[[53]],Ntimes,pg*pg);
intZHZt<-list();
for (i in 1:Ntimes) { intZHZt[[i]]<-matrix(IintZHZt[i,],pg,pg);
colnames(intZHZt[[i]])<-rownames(intZHZt[[i]])<-namesZ;  }

pscoreiid<-matrix(semiout[[52]],Ntimes,fdata$antclust*pg); 
pscore.iid<-list(); 
for (i in (0:(fdata$antclust-1))*pg) {
pscore.iid[[(i/pg)+1]]<-as.matrix(pscoreiid[,i+(1:pg)]);
colnames(pscore.iid[[i/pg+1]])<-namesZ; 
}

obs.pscore<-matrix(semiout[[51]],Ntimes,pg+1); 
sup.pscore<-as.vector(apply(as.matrix(abs(obs.pscore[,-1])),2,max)); 

ptest<-matrix(semiout[[54]],pseudo.score,pg); 

pstest.pval<-c()
for (i in 1:pg) pstest.pval<-c(pstest.pval,pval(ptest[,i],sup.pscore[i]))
}
else {sup.pscore<-intZHZt<-pscore.iid<-obs.pscore<-pstest.pval<-NULL}


if (covariance==1)  {
covit<-matrix(semiout[[35]],Ntimes,px*px);
cov.list<-list();
for (i in 1:Ntimes) cov.list[[i]]<-matrix(covit[i,],px,px);
} else cov.list<-NULL;

bhat<-matrix(semiout[[14]],nb,px+1); cum <-matrix(semiout[[15]],Ntimes,px+1); 
Vcum <-matrix(semiout[[16]],Ntimes,px+1); 
robVcum <-matrix(semiout[[17]],Ntimes,px+1); robVcum[,1] <-cum[,1]; 
gamma<-matrix(semiout[[18]],pg,1); Vargam<-matrix(semiout[[19]],pg,pg); 
robVargam<-matrix(semiout[[20]],pg,pg); 
intZHZ<-matrix(semiout[[42]],pg,pg); intZHdN<-matrix(semiout[[43]],pg,1); 
intHdN<-matrix(semiout[[55]],antpers,1); 
intHZ<-matrix(semiout[[56]],antpers,pg); 
varintHdN<-matrix(semiout[[57]],antpers,antpers); 

if (retur==1) {
cumAi<-matrix(semiout[[32]],Ntimes,fdata$antpers*1); 
cumAi<-list(time=Vcum[,1],dM=cumAi); 
#cumAI<-list(); 
#for (i in (0:(fdata$antclust-1))*p) cumAI[[i/p+1]]<-cumAi[,i+(1:p)]
} else cumAi<-NULL;

unifCI<-c();
if (sim>=1) {
test<-matrix(semiout[[23]],antsim,3*px);
testOBS<-semiout[[25]];
for (i in 1:(3*px)) testval<-c(testval,pval(test[,i],testOBS[i]))
Uit<-matrix(semiout[[29]],Ntimes,50*px); 
UIt<-list(); for (i in (0:49)*px) UIt[[i/px+1]]<-Uit[,i+(1:px)];
Ut<-matrix(semiout[[28]],Ntimes,(px+1)); 
for (i in 1:px) unifCI<-as.vector(c(unifCI,percen(test[,i],0.95)));
pval.testBeq0<-as.vector(testval[1:px]);
pval.testBeqC<-as.vector(testval[(px+1):(2*px)]);
pval.testBeqC.is<-as.vector(testval[(2*px+1):(3*px)]);
obs.testBeq0<-as.vector(testOBS[1:px]);
obs.testBeqC<-as.vector(testOBS[(px+1):(2*px)]);
obs.testBeqC.is<-as.vector(testOBS[(2*px+1):(3*px)]);
sim.testBeq0<-as.matrix(test[,1:px]);
sim.testBeqC<-as.matrix(test[,(px+1):(2*px)]);
sim.testBeqC.is<-as.matrix(test[,(2*px+1):(3*px)]);
} else {test<-NULL; unifCI<-NULL;UIt<-NULL; Ut<-NULL;
pval.testBeq0<-NULL;pval.testBeqC<-NULL; 
obs.testBeq0<-NULL;obs.testBeqC<-NULL;
sim.testBeq0<-NULL; sim.testBeqC<-NULL;
sim.testBeqC.is<-NULL; pval.testBeqC.is<-NULL; 
obs.testBeqC.is<-NULL; 
}

ud<-list(cum=cum,var.cum=Vcum,robvar.cum=robVcum,
gamma=gamma,var.gamma=Vargam,robvar.gamma=robVargam,residuals=cumAi,
pval.testBeq0=pval.testBeq0, obs.testBeq0=obs.testBeq0, 
pval.testBeqC=pval.testBeqC,pval.testBeqC.is=pval.testBeqC.is,
obs.testBeqC=obs.testBeqC, obs.testBeqC.is=obs.testBeqC.is,
sim.testBeq0= sim.testBeq0,
sim.testBeqC=sim.testBeqC,sim.testBeqC.is=sim.testBeqC.is,
conf.band=unifCI,test.procBeqC=Ut,sim.test.procBeqC=UIt,
covariance=cov.list,B.iid=B.iid,gamma.iid=gamma.iid,deviance=devi,
intZHZ=intZHZ,intZHdN=intZHdN,
obs.pscore=obs.pscore,pscore.iid=pscore.iid,intZHZt=intZHZt,
pstest.pval=pstest.pval,sup.pscore=sup.pscore,
intHdN=intHdN,intHZ=intHZ,varintHdN=varintHdN)

return(ud);
}
