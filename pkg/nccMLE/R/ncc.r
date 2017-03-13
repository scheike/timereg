
em.ncc<-function(Zobs,Tobs,status,tau,nno,Nit=1,beta=0,detail=0,
betait=2,em.dif=1,stratum=0,baseline=0,emvar=0,fullobsZ=0,
baseline.stratum=0,stratZ=0)
{
Tno<-rep(tau,nno); nobs<-length(Tobs); ntot<-nno+nobs;
ord<-order(Tobs); Tobs<-sort(Tobs); 
Zobs<-as.matrix(Zobs); Zobs[1:nobs,]<-as.matrix(Zobs)[ord,];
status<-status[ord]; px<-ncol(as.matrix(Zobs)); 
nx<-nrow(as.matrix(Zobs)); 

if (nx!=nobs) print("Du styrer for vildt"); 
if (length(tau)==nno) indi<-1 else  indi<-0; 
if (indi==0) taul<-rep(tau,nno) else taul<-tau; 

pz<-0; fullZ<-0; 
if (sum(abs(fullobsZ))!=0) {
fullZ<-as.matrix(fullobsZ); pz<-ncol(fullZ);  if (indi==0) tau<-taul; indi<-1;}

if (sum(abs(beta))==0) betaS<-rep(0,px) else betaS<-beta; loglike<-0;
loglike<-0; varbeta<-matrix(0,px,px); delta<-score<-betaS;

if (length(baseline.stratum)>1) {
baseline.stratum[1:nobs]<-baseline.stratum[ord]; antbase<-length(unique(baseline.stratum)); 
baseline.stratum<-as.numeric(factor(baseline.stratum),labels=1:antbase); 
if (length(tau)==1) tau<-rep(tau,nno); indi<-1; 
} else {antbase<-1; baseline.stratum<-rep(1,ntot);}

if (length(stratum)==ntot) {
stratum[1:nobs]<-stratum[ord]; antstrat<-length(unique(stratum)); 
stratum<-as.numeric(factor(stratum),labels=1:antstrat); 
if (length(tau)==1) tau<-rep(tau,nno); indi<-1; 
} else {antstrat<-1; nstrat<-ntot; stratum<-rep(1,ntot);}

if (stratZ==1) {z<-as.matrix(baseline.stratum[1:nobs])} else {
z<-as.matrix(Zobs); z<-as.matrix(z[,1:(px-pz)]); } 
keys<-apply(z,1,function(x) paste(x, collapse="\r"))
z<-as.matrix(z[!duplicated(keys),]); nuobs<-nrow(z);
keysz<-apply(z,1,function(x) paste(x, collapse="\r"))
orderz<-order(keysz); z<-as.matrix(z[orderz,]);
oWj<-as.matrix(table(keys,stratum[1:nobs]));
nstrat<-table(stratum); 
scases<-apply(oWj,2,sum); 

#if (win==0) dyn.load("cc.so")  else dyn.load("cc.dll");
#dyn.load("allcc.o") # LOAD NECESSARY ROUTINES

nbre<-sum(status); 
if (is.matrix(baseline)==F) breslow<-matrix(0,nbre,antbase+1) else breslow<-baseline; 
if (ncol(breslow)!=(antbase+1)) print("baseline dim ikke antal stratum"); 
breslow[,1]<-Tobs[status==1]; antW<-nrow(z);
for (i in 2:(antbase+1)) { 
breslow[,i]<-cumsum(rep(1/nbre,nbre))*max(Tobs)*
sum(status[baseline.stratum[1:nobs]==(i-1)])/(sum(c(Tobs,taul)[baseline.stratum==(i-1)])); 
} 
#print(nstrat);  print(table(stratum[1:nobs],status[1:nobs]))
#print(oWj); 

p<-matrix(0,antW,antstrat);
for (i in 1:antstrat) p[,i]<-oWj[,i]/scases[i]; #print(apply(p,2,sum)); 
konv<-rep(0,3); konv1<-rep(0,3); 
alphaij<-matrix(0,nno,antW); 

#print("bruger fast p"); 
#p<-prop.table(table(zorig,stratumo),2)
#
#const<-0.35; 
#    xx <-breslow[,1];  
#    yy <- ifelse(xx< (pi/2), const*sin(2*xx), sin(2*xx));  
#    delta <- diff(xx[1:2])
#haz<-exp(base+1*yy); 
#breslow<-cbind(xx, delta*cumsum(haz),delta*cumsum(haz*exp(1*yy)),
#delta*cumsum(haz*exp(2*yy)))
# starter med den rigtige baseline og p

if (antstrat>1 & indi==0) {indi<-1; tau<-rep(tau,nno);}

if (indi==0)  {cat("Censoring the same for all subjects that are not cases and controls \n"); 

print(c(nobs,nno,nbre,antW,ntot)); 

nparcox<-.C("emCC",
as.double(Zobs),as.integer(nobs),as.integer(px),
as.double(Tobs),as.integer(status),as.double(tau),
as.integer(nno),as.double(betaS),as.integer(Nit),
as.double(loglike),as.double(varbeta),as.double(score),
as.integer(detail),as.double(breslow),as.integer(nbre),
as.double(p),as.double(z),as.integer(antW),
as.double(oWj),as.integer(betait),as.double(delta),
as.integer(ntot),PACKAGE="nccMLE");} else {
nparcox<- .C("emCCindi",
as.double(Zobs),as.integer(nx),as.integer(px),
as.double(Tobs),as.integer(status),as.double(tau),
as.integer(nno),as.double(betaS),as.integer(Nit),
as.double(loglike),as.double(varbeta),as.double(score),
as.integer(detail),as.double(breslow),as.integer(nbre),
as.double(p),as.double(z),as.integer(antW),
as.double(oWj),as.integer(betait),as.double(delta),
as.integer(ntot),as.integer(antstrat),as.integer(stratum),
as.integer(nstrat),as.double(fullZ),as.integer(pz),
as.integer(antbase),as.integer(baseline.stratum),as.double(konv),
as.integer(stratZ),as.double(alphaij),PACKAGE="nccMLE");
konv<-nparcox[[30]]; 
nstrat<-nparcox[[25]]
alphaij<-matrix(nparcox[[32]],nno,antW); 
} 
alphany<-matrix(0,nno,antW); 

beta<-nparcox[[8]]; score<-nparcox[[12]]; 
breslow<-matrix(nparcox[[14]],nbre,antbase+1); 
p<-matrix(nparcox[[16]],antW,antstrat); delta<-nparcox[[21]];  

if (emvar==1) {
em.dif<-em.dif/ntot; 

varem<- .C("varemCC",
as.double(Zobs),as.integer(nobs),as.integer(px),as.double(Tobs),
as.integer(status),as.double(tau),as.integer(nno),as.double(beta),
as.integer(Nit),as.double(loglike),as.double(varbeta),as.double(score),
as.integer(detail),as.double(breslow),as.integer(nbre),as.double(p),
as.double(z),as.integer(antW),as.double(oWj),as.integer(betait),
as.double(delta),as.integer(ntot),as.integer(antstrat),
as.integer(stratum),as.integer(nstrat),as.double(em.dif),as.integer(indi),
as.double(fullZ),as.integer(pz),as.integer(antbase),as.integer(baseline.stratum),
as.double(konv1),as.integer(stratZ),as.double(alphany),PACKAGE="nccMLE");

d2l<--matrix(varem[[11]],px,px); 
varbeta<-solve(d2l); 
}

if (emvar==2) {
###library(numDeriv)

scoref<-function(beta,indi=0) 
{
Nit=1; betait=1; 
if (indi==0)  {
score<-.C("emCC",
as.double(Zobs),as.integer(nobs),as.integer(px),
as.double(Tobs),as.integer(status),as.double(tau),
as.integer(nno),as.double(beta),as.integer(Nit),
as.double(loglike),as.double(varbeta),as.double(score),
as.integer(detail),as.double(breslow),as.integer(nbre),
as.double(p),as.double(z),as.integer(antW),
as.double(oWj),as.integer(betait),as.double(delta),
as.integer(ntot),PACKAGE="nccMLE")[[12]];
}
else {
score <- .C("emCCindi",
as.double(Zobs),as.integer(nx),as.integer(px),
as.double(Tobs),as.integer(status),as.double(tau),
as.integer(nno),as.double(beta),as.integer(Nit),
as.double(loglike),as.double(varbeta),as.double(score),
as.integer(detail),as.double(breslow),as.integer(nbre),
as.double(p),as.double(z),as.integer(antW),
as.double(oWj),as.integer(betait),as.double(delta),
as.integer(ntot),as.integer(antstrat),as.integer(stratum),
as.integer(nstrat),as.double(fullZ),as.integer(pz),
as.integer(antbase),as.integer(baseline.stratum),as.double(konv),
as.integer(stratZ),as.double(alphaij),PACKAGE="nccMLE")[[12]];
} 
return(score)
}

###print(score)
###print(scoref(beta) )
###print(scoref(beta+c(1,0)*0.02))/0.02
###print(scoref(beta+c(0,1)*0.02))/0.02

###d2l<- -jacobian(scoref,beta,method.args=list(eps=1e-2, d=0.001)); 
###print(d2l)
###varbeta<-solve(d2l); 
###print(varbeta)
###varbeta<- d2l; 
###print(varbeta)
}

ud<-list(beta=beta,baseline=breslow,p=p,varbeta=varbeta,
delta=delta,em.dif=em.dif,covz=z,score=score,konv=konv)
return(ud)
}
