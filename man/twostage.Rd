\name{twostage}
\alias{twostage}
\title{Fits Clayton-Oakes or bivariate Plackett models for bivariate survival data
using marginals that are on Cox or addtive form.
If clusters contain more than two times, the algoritm uses a compososite likelihood
based on the pairwise bivariate models.}
\usage{
twostage(margsurv, data = sys.parent(), score.method = "nlminb", Nit = 60,
  detail = 0, clusters = NULL, silent = 1, weights = NULL,
  control = list(), theta = NULL, theta.des = NULL, var.link = 1,
  iid = 1, step = 0.5, notaylor = 0, model = "plackett",
  marginal.trunc = NULL, marginal.survival = NULL, marginal.status = NULL,
  strata = NULL, se.clusters = NULL, max.clust = NULL, numDeriv = 1)
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

  \item{marginal.trunc}{marginal left truncation
  probabilities}

  \item{marginal.survival}{optional vector of marginal
  survival probabilities}

  \item{marginal.status}{related to marginal survival
  probabilities}

  \item{strata}{strata for fitting, see example}

  \item{se.clusters}{for clusters for se calculation with
  iid}

  \item{max.clust}{max se.clusters for se calculation with
  iid}

  \item{numDeriv}{to get numDeriv version of second
  derivative, otherwise uses sum of squared score}
}
\description{
The reported standard errors are based on the estimated
information from the likelihood assuming that the marginals
are known.
}
\examples{
data(diabetes)

# Marginal Cox model  with treat as covariate
margph <- coxph(Surv(time,status)~treat,data=diabetes)
### Clayton-Oakes, from timereg
fitco1<-two.stage(margph,data=diabetes,theta=1.0,detail=0,Nit=40,clusters=diabetes$id)
summary(fitco1)
### Plackett model
fitp<-twostage(margph,data=diabetes,theta=0.0,Nit=40,
               clusters=diabetes$id,var.link=1)
summary(fitp)
### Clayton-Oakes
fitco2<-twostage(margph,data=diabetes,theta=1.0,detail=0,Nit=40,
                 clusters=diabetes$id,var.link=1,model="clayton.oakes")
summary(fitco2)

### without covariates using Aalen for marginals
marg <- aalen(Surv(time,status)~+1,data=diabetes,n.sim=0,max.clust=NULL,robust=0)
fitpa<-twostage(marg,data=diabetes,theta=1.0,detail=0,Nit=40,
                clusters=diabetes$id,score.method="optimize")
summary(fitpa)

fitcoa<-twostage(marg,data=diabetes,theta=1.0,detail=0,Nit=40,clusters=diabetes$id,
                 var.link=1,model="clayton.oakes")
summary(fitcoa)

### Piecewise constant cross hazards ratio modelling
########################################################

d <- subset(simClaytonOakes(2000,2,0.5,0,stoptime=2,left=0),!truncated)
udp <- piecewise.twostage(c(0,0.5,2),data=d,score.method="optimize",
                          id="cluster",timevar="time",
                          status="status",model="clayton.oakes",silent=0)
summary(udp)

\donttest{
### Same model using the strata option, a bit slower
########################################################
## makes the survival pieces for different areas in the plane
##ud1=surv.boxarea(c(0,0),c(0.5,0.5),data=d,id="cluster",timevar="time",status="status")
##ud2=surv.boxarea(c(0,0.5),c(0.5,2),data=d,id="cluster",timevar="time",status="status")
##ud3=surv.boxarea(c(0.5,0),c(2,0.5),data=d,id="cluster",timevar="time",status="status")
##ud4=surv.boxarea(c(0.5,0.5),c(2,2),data=d,id="cluster",timevar="time",status="status")

## everything done in one call
ud <- piecewise.data(c(0,0.5,2),data=d,timevar="time",status="status",id="cluster")
ud$strata <- factor(ud$strata);
ud$intstrata <- factor(ud$intstrata)

## makes strata specific id variable to identify pairs within strata
## se's computed based on the id variable across strata "cluster"
ud$idstrata <- ud$id+(as.numeric(ud$strata)-1)*2000

marg2 <- aalen(Surv(boxtime,status)~-1+factor(num):factor(intstrata),
               data=ud,n.sim=0,robust=0)
tdes <- model.matrix(~-1+factor(strata),data=ud)
fitp2<-twostage(marg2,data=ud,se.clusters=ud$cluster,clusters=ud$idstrata,
                score.method="fisher.scoring",model="clayton.oakes",
                theta.des=tdes,step=0.5)
summary(fitp2)

### now fitting the model with symmetry, i.e. strata 2 and 3 same effect
ud$stratas <- ud$strata;
ud$stratas[ud$strata=="0.5-2,0-0.5"] <- "0-0.5,0.5-2"
tdes2 <- model.matrix(~-1+factor(stratas),data=ud)
fitp3<-twostage(marg2,data=ud,clusters=ud$idstrata,se.cluster=ud$cluster,
                score.method="fisher.scoring",model="clayton.oakes",
                theta.des=tdes2,step=0.5)
summary(fitp3)

### same model using strata option, a bit slower
fitp4<-twostage(marg2,data=ud,clusters=ud$cluster,se.cluster=ud$cluster,
                score.method="fisher.scoring",model="clayton.oakes",
                theta.des=tdes2,step=0.5,strata=ud$strata)
summary(fitp4)
}
}
\author{
Thomas Scheike
}
\references{
Clayton-Oakes and Plackett bivariate survival
distributions,
}
\keyword{survival}

