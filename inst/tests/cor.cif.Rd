library(mets)
library(vimcom)
data(multcif);
multcif$cause[multcif$cause==0] <- 2
zyg <- rep(rbinom(200,1,0.5),each=2)
theta.des <- model.matrix(~-1+factor(zyg))

times=seq(0.05,1,by=0.05) # to speed up computations use only these time-points
add<-comp.risk(Surv(time,cause>0)~+1+cluster(id),data=multcif,causeS=1,
               multcif$cause,n.sim=0,times=times,model="fg",max.clust=NULL)
add2<-comp.risk(Surv(time,cause>0)~+1+cluster(id),data=multcif,causeS=2,
               multcif$cause,n.sim=0,times=times,model="fg",max.clust=NULL)

out1<-cor.cif(add,data=multcif,cause1=1,cause2=1)
summary(out1)

out2<-cor.cif(add,data=multcif,cause1=1,cause2=1,theta.des=theta.des)
summary(out2)

out3<-cor.cif(add,data=multcif,cause1=1,cause2=2,cif2=add2)
summary(out3)
###########################################################
# investigating further models using parfunc and dparfunc
###########################################################
set.seed(100)
prt<-simnordic(1000,cordz=2,cormz=5)
prt$status <-prt$cause
table(prt$status)

times <- seq(40,100,by=10)
cifmod <- comp.risk(Surv(time,status>0)~+1+cluster(id),data=prt,prt$status,causeS=1,n.sim=0,
                  times=times,conservative=1,max.clust=NULL,model="fg")
theta.des <- model.matrix(~-1+factor(zyg),data=prt)

parfunc <- function(par,t,pardes)
{
par <- pardes %*% c(par[1],par[2]) +
       pardes %*% c( par[3]*(t-60)/12,par[4]*(t-60)/12)
par
}
head(parfunc(c(0.1,1,0.1,1),50,theta.des))

dparfunc <- function(par,t,pardes)
{
dpar <- cbind(pardes, t(t(pardes) * c( (t-60)/12,(t-60)/12)) )
dpar
}
head(dparfunc(c(0.1,1,0.1,1),50,theta.des))

names(prt)
or1 <- or.cif(cifmod,data=prt,cause1=1,cause2=1,theta.des=theta.des,same.cens=TRUE,theta=c(0.6,1.1,0.1,0.1),
  par.func=parfunc,dpar.func=dparfunc,dimpar=4,score.method="fisher.scoring",detail=1)
summary(or1)

 cor1 <- cor.cif(cifmod,data=prt,cause1=1,cause2=1,theta.des=theta.des,same.cens=TRUE,theta=c(0.5,1.0,0.1,0.1),
       par.func=parfunc,dpar.func=dparfunc,dimpar=4,control=list(trace=TRUE),detail=1)
summary(cor1)

### piecwise contant OR model
gparfunc <- function(par,t,pardes)
{
	cuts <- c(0,80,90,120)
	grop <- diff(t<cuts)
paru  <- (pardes[,1]==1) * sum(grop*par[1:3])+ (pardes[,2]==1) * sum(grop*par[4:6])
paru
}

dgparfunc <- function(par,t,pardes)
{
	cuts <- c(0,80,90,120)
	grop <- diff(t<cuts)
par1 <- matrix(c(grop),nrow(pardes),length(grop),byrow=TRUE)
parmz <- par1* (pardes[,1]==1)
pardz <- (pardes[,2]==1) * par1
dpar <- cbind( parmz,pardz)
dpar
}
head(dgparfunc(rep(0.1,6),50,theta.des))
head(gparfunc(rep(0.1,6),50,theta.des))

or1g <- or.cif(cifmod,data=prt,cause1=1,cause2=1,theta.des=theta.des,same.cens=TRUE,
       par.func=gparfunc,dpar.func=dgparfunc,dimpar=6,score.method="fisher.scoring",detail=1)
summary(or1g)
names(or1g)
or1g$theta.iid

