
library(mets)

set.seed(100)
n <- 40000

## {{{ competing risks ace model with profile of baseline 
lam0 <- c(0.3,0.2)
pars <- c(1,1,1,1,0.1,1)*0.25
## genetic random effects, cause1, cause2 and overall 
parg <- pars[c(1,3,5)]
## environmental random effects, cause1, cause2 and overall 
parc <- pars[c(2,4,6)]

## simulate competing risks with two causes with hazards 0.5 and 0.3
## ace for each cause, and overall ace 
out <- simCompete.twin.ace(n,parg,parc,0,2,
			   lam0=lam0,overall=1,all.sum=1)

## setting up design for running the model 
## {{{ setting pairs and random effects 
# 
mm <- familycluster.index(out$cluster)
head(mm$familypairindex,n=10)
pairs <- matrix(mm$familypairindex,ncol=2,byrow=TRUE)
tail(pairs,n=12)
#
kinship <- (out[pairs[,1],"zyg"]=="MZ")+ (out[pairs[,1],"zyg"]=="DZ")*0.5

dout <- make.pairwise.design.competing(pairs,kinship,
	       type="ace",compete=length(lam0),overall=1)
head(dout$ant.rvs)
## MZ
dim(dout$theta.des)
dout$theta.des[,,1]
dout$random.design[,,1]
dout$theta.des[,,nrow(out)/2]
dout$random.design[,,nrow(out)/2]
###table(out$status)
## }}} 

## competing risks models, given as list 
cr.models=list(Surv(time,status==1)~+1,
	       Surv(time,status==2)~+1)
ms <- out$time %o% lam0

par(mfrow=c(1,1))
tsf <- twostage(NULL,data=out,clusters=out$cluster,
               theta=0.01+pars/sum(pars)^2,
	       var.link=0,step=1.0,Nit=20,detail=1,
               random.design=dout$random.design,
               theta.des=dout$theta.des,pairs=pairs,
	       numDeriv=0,
               marginal.status=out$status,
               marginal.survival=ms,
	       two.stage=0,cr.models=cr.models)
coef(tsf)
pars/sum(pars)^2
summary(tsf$marginal.trunc)
summary(tsf$marginal.surv)
###tsf$score; tsf$score1

###
source("../R/twostage.R")
par(mfrow=c(1,1))
ts <- twostage(NULL,data=out,clusters=out$cluster,
               theta=tsf$theta,
	       step=1.0,Nit=20,detail=1,
               random.design=dout$random.design,
               theta.des=dout$theta.des,pairs=pairs,
	       numDeriv=0,
               marginal.status=out$status,
	       two.stage=0,cr.models=cr.models)

coef(ts)
pars/sum(pars)^2

matplot.twostage(ts)
abline(c(0,lam0[1])); abline(c(0,lam0[2])); 

system.time(
out1 <- aalen(cr.models[[1]],data=out,robust=0)
)
system.time(
out1 <- aalen(Surv(time,status!=0)~+1,data=out,robust=0)
)


## }}} 

onec <- 1
if (onec==1) {
## {{{ one cause ACE survival
# 

lam0 <- c(0.5)
pars <- c(1,0.5); 
pars <- c(0.5,1); 
out <- simCompete.twin.ace(n,pars[1],pars[2],0,2,lam0=lam0,overall=0)

## {{{ 
mm <- familycluster.index(out$cluster)
head(mm$familypairindex,n=10)
pairs <- matrix(mm$familypairindex,ncol=2,byrow=TRUE)
tail(pairs,n=12)
#
kinship <- (out[pairs[,1],"zyg"]=="MZ")+ (out[pairs[,1],"zyg"]=="DZ")*0.5

dout <- make.pairwise.design(pairs,kinship,type="ace")
head(dout$ant.rvs)
## MZ
dim(dout$theta.des)
dout$theta.des[,,1]
dout$random.design[,,1]
## DZ
dout$theta.des[,,nrow(out)/2]
dout$random.design[,,nrow(out)/2]
## }}} 

lams <- cbind(lam0[1]*out$time)
table(out$status)

ts <- twostage(NULL,data=out,clusters=out$cluster,
               theta=pars/sum(pars)^2,
	       var.link=0,step=1.0,Nit=10,detail=0,
               random.design=dout$random.design,
               theta.des=dout$theta.des,pairs=pairs,
	       numDeriv=0,
               marginal.surv=lams,
	       marginal.status=out$status,
               two.stage=0)
summary(ts)

ts2 <- twostage(NULL,data=out,clusters=out$cluster,
               theta=ts$theta,
	       var.link=0,step=1.0,Nit=10,detail=1,
               random.design=dout$random.design,
               theta.des=dout$theta.des,pairs=pairs,
	       numDeriv=0,
	       marginal.status=out$status,
	       cr.model=list(Surv(time,status)~+1),
               two.stage=0)
summary(ts2)
ts2$score

## }}} 
}

twoc <- 1
if (twoc==1) {
## {{{ two cause two independent ACE survival
# 

lam0 <- c(0.5,0.4)
pars <- c(0.5,1,0.5,1)*0.5; 
out <- simCompete.twin.ace(n,pars[c(1,3)],pars[c(2,4)],0,2,lam0=lam0,overall=0)
table(out$status)
out$status1 <- out$status==1
table(out$status1)

## {{{ 
mm <- familycluster.index(out$cluster)
head(mm$familypairindex,n=10)
pairs <- matrix(mm$familypairindex,ncol=2,byrow=TRUE)
tail(pairs,n=12)
#
kinship <- (out[pairs[,1],"zyg"]=="MZ")+ (out[pairs[,1],"zyg"]=="DZ")*0.5

###kinship <- c(rep(1,5000),rep(0.5,5000))
###
###source("mets/R/twostage.R")
###dout <- make.pairwise.design(pairs,kinship,type="ace")
dout <- make.pairwise.design.competing(pairs,kinship,
	       type="ace",overall=0,compete=2)
head(dout$ant.rvs)
## MZ
dim(dout$theta.des)
dout$theta.des[,,1]
dout$random.design[,,1]
dout$theta.des[,,nrow(out)/2]
dout$random.design[,,nrow(out)/2]
###
###
###out$status[out$status==3] <- 2
# 
table(out$status)
## }}} 

###lams <- cbind(lam0[1]*out$time,lam0[2]*out$time)
lams <- cbind(out$time*lam0[1],out$time*lam0[2])

ts <- twostage(NULL,data=out,clusters=out$cluster,
               theta=pars/sum(pars)^2,
	       var.link=0,step=1.0,Nit=10,detail=1,
               random.design=dout$random.design,
               theta.des=dout$theta.des,pairs=pairs,
               marginal.surv=lams,
	       marginal.status=out$status,
               two.stage=0)
summary(ts)
ts$score

ts2 <- twostage(NULL,data=out,clusters=out$cluster,
               theta=pars/sum(pars)^2,
	       var.link=0,step=1.0,Nit=20,detail=1,
               random.design=dout$random.design,
               theta.des=dout$theta.des,pairs=pairs,
	       marginal.status=out$status,
               two.stage=0,
     cr.models=list( Surv(time,status==1)~+1,
		     Surv(time,status==2)~+1)
	       )
ts2$theta
pars/sum(pars)^2

matplot.twostage(ts2)
abline(c(0,lam0[1])); abline(c(0,lam0[2])); 

} ## }}} 


itwoc <- 1
if (itwoc==1) {
## {{{ cause specific analyses because independence 

lam0 <- c(0.5,0.4)
pars <- c(1,0.5); 
pars <- c(0.5,1,0.5,1); 
out <- simCompete.twin.ace(n,
	   pars[c(1,3)],pars[c(2,4)],
           0,2,lam0=lam0,overall=0,all.sum=1)
table(out$status)

## {{{ 
mm <- familycluster.index(out$cluster)
head(mm$familypairindex,n=10)
pairs <- matrix(mm$familypairindex,ncol=2,byrow=TRUE)
tail(pairs,n=12)
#
kinship <- (out[pairs[,1],"zyg"]=="MZ")+ (out[pairs[,1],"zyg"]=="DZ")*0.5

###kinship <- c(rep(1,5000),rep(0.5,5000))
###
###source("mets/R/twostage.R")
###dout <- make.pairwise.design(pairs,kinship,type="ace")
dout <- make.pairwise.design.competing(pairs,kinship,
	       type="ace",overall=0)
head(dout$ant.rvs)
## MZ
dim(dout$theta.des)
dout$theta.des[,,1]
dout$random.design[,,1]
## DZ
dout$theta.des[,,nrow(out)/2]
dout$random.design[,,nrow(out)/2]
###
###
# 
table(out$status)
out$statusc1 <- 1*(out$status==1)
out$statusc2 <- 1*(out$status==2)
## }}} 


## design for only cause 1
dout1 <- make.pairwise.design(pairs,kinship,
	       type="ace")

## competing risks models, given as list 
cr.models=list(Surv(time,status==1)~+1,Surv(time,status==2)~+1)

ts <- twostage(NULL,data=out,clusters=out$cluster,
               theta=pars,
	       score.method="fisher.scoring",
	       var.link=0,step=1.0,Nit=10,detail=0,
               random.design=dout$random.design,
               theta.des=dout$theta.des,pairs=pairs,
	       numDeriv=0,
               marginal.status=out$status,
       two.stage=0,cr.models=cr.models)
coef(ts)


### considering cause 1 alone 

cr.models1=list(Surv(time,statusc1==1)~+1)

## note due to parametrization ags=2 ! 
ts1 <- twostage(NULL,data=out,clusters=out$cluster,
               theta=pars[1:2],
	       var.link=0,step=1.0,Nit=10,detail=0,
               random.design=dout1$random.design,
               theta.des=dout1$theta.des,
	       pairs=pairs,
	       marginal.status=out$statusc1,
               two.stage=0,
	       cr.models=cr.models2)
summary(ts1)

### considering cause 2 alone 

cr.models2=list(Surv(time,statusc2==1)~+1)

ts2 <- twostage(NULL,data=out,clusters=out$cluster,
               theta=pars[1:2],
	       var.link=0,step=1.0,Nit=10,detail=0,
               random.design=dout1$random.design,
               theta.des=dout1$theta.des,pairs=pairs,
	       marginal.status=out$statusc2,
               two.stage=0,
	       cr.models=cr.models2)
summary(ts2)

## }}} 
}


## {{{ 2 compete+ overall ace, dependence test, cox 
### conditional cox model 
lam0 <- c(0.5,0.5)
pars <- rep(1,6)
out <- simCompete.twin.ace(n,c(1,1,1),c(1,1,1),0,2,
		  lam0=lam0,overall=1,all.sum=1)

table(out$status)
###
out2 <- fast.reshape(out,id="cluster")
###
out2$mintime <- pmin(out2$time1,out2$time2)
out2$whichmin <- ifelse(out2$time1<out2$time2,1,2)
out2$mincause <- ifelse(out2$time1<out2$time2,out2$status1,
			out2$status2)
out2$maxtime <- pmax(out2$time1,out2$time2)
out2$maxcause <- ifelse(out2$time1>out2$time2,out2$status1,
			out2$status2)
###
out1 <- event.split(out2,time="time1",status="status1",
	    cuts="time2")
out1$mstat <- out1$status2*out1$num
head(out1[out1$whichmin==2,])
head(out1[out1$whichmin==1,])
table(out1$status1)
table(out1$mstat)
###out1[out1$num==1,]
out1$mincause[out1$num==0] <- 0
###
coxph(Surv(start,time1,status1==1)~factor(mstat)*factor(zyg1),data=out1)
###
coxph(Surv(start,time1,status1==2)~factor(mstat)*factor(zyg1),data=out1)

mm <- familycluster.index(out$cluster)
head(mm$familypairindex,n=10)
pairs <- matrix(mm$familypairindex,ncol=2,byrow=TRUE)
tail(pairs,n=12)
#
kinship <- (out[pairs[,1],"zyg"]=="MZ")+ (out[pairs[,1],"zyg"]=="DZ")*0.5

dout <- make.pairwise.design.competing(pairs,kinship,
	       type="ace",compete=length(lam0),overall=1)
head(dout$ant.rvs)
## MZ
dim(dout$theta.des)
dout$theta.des[,,1]
dout$random.design[,,1]
## DZ
dout$theta.des[,,nrow(out)/2]
dout$random.design[,,nrow(out)/2]
###
###

cr.models=list(Surv(time,status==1)~+1,Surv(time,status==2)~+1)

ts <- twostage(NULL,data=out,clusters=out$cluster,
               theta=pars,
	       score.method="fisher.scoring",
	       var.link=0,step=1.0,Nit=20,detail=0,
               random.design=dout$random.design,
               theta.des=dout$theta.des,pairs=pairs,
               marginal.status=out$status,
	       two.stage=0,cr.models=cr.models)

## }}} 



