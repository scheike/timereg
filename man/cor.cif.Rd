\name{cor.cif}
\alias{cor.cif}
\alias{or.cif}
\alias{rr.cif}
\title{Cross-odds-ratio, OR or RR risk regression for competing risks}
\usage{
  cor.cif(cif, data, cause, times = NULL, cause1 = 1,
    cause2 = 1, cens.code = 0, cens.model = "KM", Nit = 40,
    detail = 0, clusters = NULL, theta = NULL,
    theta.des = NULL, step = 1, sym = 0, weights = NULL,
    par.func = NULL, dpar.func = NULL, dimpar = NULL,
    score.method = "nlminb", same.cens = FALSE,
    censoring.probs = NULL, silent = 1, ...)
}
\arguments{
  \item{cif}{a model object from the comp.risk function
  with the marginal cumulative incidence of cause1, i.e.,
  the event of interest, and whose odds the comparision is
  compared to the conditional odds given cause2}

  \item{data}{a data.frame with the variables.}

  \item{cause}{specifies the causes related to the death
  times, the value cens.code is the censoring value. When
  missing it comes from marginal cif.}

  \item{times}{time-vector that specifies the times used
  for the estimating euqations for the cross-odds-ratio
  estimation.}

  \item{cause1}{specificies the cause considered.}

  \item{cause2}{specificies the cause that is conditioned
  on.}

  \item{cens.code}{specificies the code for the censoring.}

  \item{cens.model}{specified which model to use for the
  ICPW, KM is Kaplan-Meier alternatively it may be "cox"}

  \item{Nit}{number of iterations for Newton-Raphson
  algorithm.}

  \item{detail}{if 0 no details are printed during
  iterations, if 1 details are given.}

  \item{clusters}{specifies the cluster structure.}

  \item{theta}{specifies starting values for the
  cross-odds-ratio parameters of the model.}

  \item{theta.des}{specifies a regression design for the
  cross-odds-ratio parameters.}

  \item{step}{specifies the step size for the
  Newton-Raphson algorithm.}

  \item{sym}{specifies if symmetry is used in the model.}

  \item{weights}{weights for estimating equations.}

  \item{par.func}{parfunc}

  \item{dpar.func}{dparfunc}

  \item{dimpar}{dimpar}

  \item{score.method}{"nlminb", can also use
  "fisher-scoring".}

  \item{same.cens}{if true then censoring within clusters
  are assumed to be the same variable, default is
  independent censoring.}

  \item{censoring.probs}{if cens.model is "user.weights"
  these probabilities are used for the bivariate censoring
  dist.}

  \item{silent}{1 to suppress output about convergence
  related issues.}

  \item{...}{Not used.}
}
\value{
  returns an object of type 'cor'. With the following
  arguments: \item{theta}{estimate of proportional odds
  parameters of model.} \item{var.theta}{variance for
  gamma.  } \item{hess}{the derivative of the used score.}
  \item{score}{scores at final stage.} \item{score}{scores
  at final stage.} \item{theta.iid}{matrix of iid
  decomposition of parametric effects.}
}
\description{
  Fits a parametric model for the log-cross-odds-ratio for
  the predictive effect of for the cumulative incidence
  curves for \eqn{T_1} experiencing cause i given that
  \eqn{T_2} has experienced a cause k : \deqn{
  \log(COR(i|k)) = h(\theta,z_1,i,z_2,k,t)=_{default}
  \theta^T z = } with the log cross odds ratio being \deqn{
  COR(i|k) = \frac{O(T_1 \leq t,cause_1=i | T_2 \leq
  t,cause_2=k)}{ O(T_1 \leq t,cause_1=i)} } the conditional
  odds divided by the unconditional odds, with the odds
  being, respectively \deqn{ O(T_1 \leq t,cause_1=i | T_2
  \leq t,cause_1=k) = \frac{ P_x(T_1 \leq t,cause_1=i | T_2
  \leq t,cause_2=k)}{ P_x((T_1 \leq t,cause_1=i)^c | T_2
  \leq t,cause_2=k)} } and \deqn{ O(T_1 \leq t,cause_1=i) =
  \frac{P_x(T_1 \leq t,cause_1=i )}{P_x((T_1 \leq
  t,cause_1=i)^c )}. } Here \eqn{B^c} is the complement
  event of \eqn{B}, \eqn{P_x} is the distribution given
  covariates (\eqn{x} are subject specific and \eqn{z} are
  cluster specific covariates), and \eqn{h()} is a function
  that is the simple identity \eqn{\theta^T z} by default.
}
\details{
  The OR dependence measure is given by \deqn{ OR(i,k) =
  \log ( \frac{O(T_1 \leq t,cause_1=i | T_2 \leq
  t,cause_2=k)}{ O(T_1 \leq t,cause_1=i) | T_2 \leq
  t,cause_2=k)} } This measure is numerically more stabile
  than the COR measure, and is symetric in i,k.

  The RR dependence measure is given by \deqn{ RR(i,k) =
  \log ( \frac{P(T_1 \leq t,cause_1=i , T_2 \leq
  t,cause_2=k)}{ P(T_1 \leq t,cause_1=i) P(T_2 \leq
  t,cause_2=k)} } This measure is numerically more stabile
  than the COR measure, and is symetric in i,k.

  The model is fitted under symmetry (sym=1), i.e., such
  that it is assumed that \eqn{T_1} and \eqn{T_2} can be
  interchanged and leads to the same cross-odd-ratio (i.e.
  \eqn{COR(i|k) = COR(k|i))}, as would be expected for
  twins or without symmetry as might be the case with
  mothers and daughters (sym=0).

  \eqn{h()} may be specified as an R-function of the
  parameters, see example below, but the default is that it
  is simply \eqn{\theta^T z}.
}
\examples{
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

times <- seq(40,100,by=2)
cifmod <- comp.risk(Surv(time,status>0)~+1+cluster(id),data=prt,prt$status,causeS=1,n.sim=0,
                  times=times,conservative=1,max.clust=NULL,model="fg")
theta.des <- model.matrix(~-1+factor(zyg),data=prt)

parfunc <- function(par,t,pardes)
{
par <- pardes \%*\% c(par[1],par[2]) +
       pardes \%*\% c( par[3]*(t-60)/12,par[4]*(t-60)/12)
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
\dontrun{
or1 <- or.cif(cifmod,data=prt,cause1=1,cause2=1,theta.des=theta.des,same.cens=TRUE,theta=c(0.6,1.1,0.1,0.1),
       par.func=parfunc,dpar.func=dparfunc,dimpar=4,score.method="fisher.scoring",detail=1)
summary(or1)

 cor1 <- cor.cif(cifmod,data=prt,cause1=1,cause2=1,theta.des=theta.des,same.cens=TRUE,theta=c(0.5,1.0,0.1,0.1),step=0.5,
       par.func=parfunc,dpar.func=dparfunc,dimpar=4,control=list(trace=TRUE),detail=1)
summary(cor1)
### piecwise contant OR model
gparfunc <- function(par,t,pardes)
{
	cuts <- c(0,70,80,90,120)
	grop <- diff(t<cuts)
par  <- (pardes[,1]==1) * sum(grop*par[1:4])+ (pardes[,2]==1) * sum(grop*par[5:8])
par
}

dgparfunc <- function(par,t,pardes)
{
	cuts <- c(0,70,80,90,120)
	grop <- diff(t<cuts)
par1 <- matrix(c(grop),nrow(pardes),4,byrow=TRUE)
parmz <- par1[,1:4]* (pardes[,1]==1)
pardz <- (pardes[,2]==1) * par1
dpar <- cbind( parmz,pardz)
dpar
}
head(dgparfunc(c(0.1,1,0.1,1),50,theta.des))
or1g <- or.cif(cifmod,data=prt,cause1=1,cause2=1,theta.des=theta.des,same.cens=TRUE,
       par.func=gparfunc,dpar.func=dgparfunc,dimpar=8,score.method="fisher.scoring",detail=1)
summary(or1g)
}
}
\author{
  Thomas Scheike
}
\references{
  Cross odds ratio Modelling of dependence for Multivariate
  Competing Risks Data, Scheike and Sun (2010), work in
  progress.

  A Semiparametric Random Effects Model for Multivariate
  Competing Risks Data, Scheike, Zhang, Sun, Jensen (2010),
  Biometrika.
}
\keyword{survival}

