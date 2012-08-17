\name{twostage}
\alias{twostage}
\title{Fits Clayton-Oakes or bivariate Plackett models for bivariate survival data
using marginals that are on Cox or addtive form.
If clusters contain more than two times, the algoritm uses a compososite likelihood
based on the pairwise bivariate models.}
\usage{
  twostage(margsurv, data = sys.parent(),
    score.method = "nlminb", Nit = 60, detail = 0,
    start.time = 0, max.time = NULL, id = NULL,
    clusters = NULL, silent = 1, weights = NULL,
    control = list(), robust = 0, theta = NULL,
    theta.des = NULL, var.link = 1, iid = 0, step = 0.5,
    notaylor = 0, model = "plackett")
}
\description{
  The reported standard errors are based on the estimated
  information from the likelihood assuming that the
  marginals are known.
}
\examples{
data(diabetes)

data(diabetes)
# Marginal Cox model  with treat as covariate
margph <- coxph(Surv(time,status)~treat,data=diabetes)
### Clayton-Oakes, from timereg
fitco1<-two.stage(margph,data=diabetes,theta=1.0,detail=0,Nit=40,clusters=diabetes$id)
summary(fitco1)
### Plackett model
fitp<-twostage(margph,data=diabetes,theta=1.0,detail=0,Nit=40,clusters=diabetes$id,var.link=1)
summary(fitp)
### Clayton-Oakes
fitco2<-twostage(margph,data=diabetes,theta=1.0,detail=0,Nit=40,clusters=diabetes$id,var.link=1,model="clayton.oakes")
summary(fitco2)

### without covariates using Aalen for marginals
marg <- aalen(Surv(time,status)~+1,data=diabetes,n.sim=0,max.clust=NULL,robust=0)
fitpa<-twostage(marg,data=diabetes,theta=1.0,detail=0,Nit=40,clusters=diabetes$id,score.method="optimize")
summary(fitpa)

fitcoa<-twostage(marg,data=diabetes,theta=1.0,detail=0,Nit=40,clusters=diabetes$id,var.link=1,model="clayton.oakes")
summary(fitcoa)

diab1 <- surv.boxarea(c(10,10),c(40,40),diabetes) ### data truncated at (10,10) and right-censored (40,40)
diab2 <- surv.boxarea(c(40,40),c(90,90),diabetes) ### data truncated at (40,40) and right-censored (90,90)

### dependence for early events
marg1 <- aalen(Surv(left,time,status)~+1,data=diab1,n.sim=0,max.clust=NULL,robust=0)
fit1<-twostage(marg1,data=diab1,clusters=diab1$id,model="clayton.oakes")
summary(fit1)

### dependence for later events
marg2 <- aalen(Surv(left,time,status)~+1,data=diab2,n.sim=0,max.clust=NULL,robust=0)
fit2<-twostage(marg2,data=diab2,clusters=diab2$id,model="clayton.oakes")
summary(fit2)
}
\author{
  Thomas Scheike
}
\references{
  Clayton-Oakes and Plackett bivariate survival
  distributions,
}
\keyword{survival}

