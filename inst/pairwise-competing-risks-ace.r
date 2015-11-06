
library(mets)
## {{{ competing risks ace model with profile of baseline 
lam0 <- c(0.5,0.3)
pars <- c(1,1,1,1,0,1)
## genetic random effects, cause1, cause2 and overall 
parg <- pars[c(1,3,5)]
## environmental random effects, cause1, cause2 and overall 
parc <- pars[c(2,4,6)]

## simulate competing risks with two causes with hazards 0.5 and 0.3
## ace for each cause, and overall ace 
out <- simCompete.twin.ace(40000,parg,parc,0,2,lam0=lam0,overall=1,all.sum=1)

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
###dim(dout$theta.des)
###dout$theta.des[,,1]
###dout$random.design[,,1]
###dout$theta.des[,,nrow(out)/2]
###dout$random.design[,,nrow(out)/2]
###
###
out$status[out$status==3] <- 2
# 
###table(out$status)
## }}} 

## competing risks models, given as list 
cr.models=list(Surv(time,status==1)~+1,Surv(time,status==2)~+1)

ts <- twostage(NULL,data=out,clusters=out$cluster,
               theta=pars,
	       score.method="fisher.scoring",
	       var.link=0,step=1.0,Nit=5,detail=0,
               random.design=dout$random.design,
               theta.des=dout$theta.des,pairs=pairs,
	       numDeriv=0,
               marginal.status=out$status,
	       two.stage=0, cr.models=cr.models)

system.time(
out1 <- aalen(cr.models[[1]],data=out,robust=0)
)
system.time(
out1 <- aalen(Surv(time,status!=0)~+1,data=out,robust=0)
)

## }}} 

