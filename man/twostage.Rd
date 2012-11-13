\name{twostage}
\alias{twostage}
\title{Fits Clayton-Oakes or bivariate Plackett models for bivariate survival data
using marginals that are on Cox or addtive form.
If clusters contain more than two times, the algoritm uses a compososite likelihood
based on the pairwise bivariate models.}
\usage{
  twostage(margsurv, data = sys.parent(),
    score.method = "nlminb", Nit = 60, detail = 0,
    clusters = NULL, silent = 1, weights = NULL,
    control = list(), theta = NULL, theta.des = NULL,
    var.link = 1, iid = 0, step = 0.5, notaylor = 0,
    model = "plackett", marginal.survival = NULL,
    strata = NULL)
}
\arguments{
  \item{margsurv}{Marginal model}

  \item{data}{data frame}

  \item{score.method}{Scoring method}

  \item{Nit}{Number of iterations}

  \item{detail}{Detail}

  \item{clusters}{Cluster variable}

  \item{silent}{Debug information}

  \item{weights}{Weights}

  \item{control}{Optimization arguments}

  \item{theta}{Starting values for variance components}

  \item{theta.des}{Variance component design}

  \item{var.link}{Link function for variance}

  \item{iid}{Calculate i.i.d. decomposition}

  \item{step}{Step size}

  \item{notaylor}{Taylor expansion}

  \item{model}{model}

  \item{marginal.survival}{optional vector of marginal
  survival probabilities}

  \item{strata}{Strata}
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

### Piecewise constant cross hazards ratio modelling
########################################################

d <- subset(simClaytonOakes(2000,2,0.5,0,stoptime=2,left=0),!truncated)
udp <- piecewise.twostage(c(0,0.5,2),data=d,score.method="optimize",id="cluster",timevar="time",
			  status="status",model="clayton.oakes",silent=0)
summary(udp)

\donttest{
### Same model using the strata option, a bit slower
########################################################

ud1=surv.boxarea(c(0,0),c(0.5,0.5),data=d,id="cluster",timevar="time",status="status")
ud2=surv.boxarea(c(0,0.5),c(0.5,2),data=d,id="cluster",timevar="time",status="status")
ud3=surv.boxarea(c(0.5,0),c(2,0.5),data=d,id="cluster",timevar="time",status="status")
ud4=surv.boxarea(c(0.5,0.5),c(2,2),data=d,id="cluster",timevar="time",status="status")
ud1$strata <- 1; ud2$strata <- 2; ud3$strata <- 3; ud4$strata <- 4
ud <- rbind(ud1,ud2,ud3,ud4)

marg2 <- aalen(Surv(boxtime,status)~-1+factor(num):factor(strata),data=ud,n.sim=0,robust=0)
tdes <- model.matrix(~-1+factor(strata),data=ud)
fitp2<-twostage(marg2,data=ud,clusters=ud$id,score.method="fisher.scoring",model="clayton.oakes",
                theta.des=tdes,step=1.0,detail=0,strata=ud$strata)
summary(fitp2)

### now fitting the model with symmetry, i.e. strata 2 and 3 same effect
ud$stratas <- ud$strata; ud$stratas[ud$strata==3] <- 2;
tdes2 <- model.matrix(~-1+factor(stratas),data=ud)
fitp3<-twostage(marg2,data=ud,clusters=ud$id,score.method="fisher.scoring",model="clayton.oakes",
                theta.des=tdes2,step=1.0,detail=0,strata=ud$strata)
summary(fitp3)
}
### could also symmetry in marginal models
}
\author{
  Thomas Scheike
}
\references{
  Clayton-Oakes and Plackett bivariate survival
  distributions,
}
\keyword{survival}

