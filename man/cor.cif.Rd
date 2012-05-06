\name{cor.cif}
\alias{cor.cif}
\alias{or.cif}
\alias{rr.cif}
\title{Cross-odds-ratio, OR or RR risk regression for competing risks}
\usage{
  cor.cif(cif, data, cause, times = NULL, cause1 = 1,
    cause2 = 1, cens.code = 0, cens.model = "KM", Nit = 40,
    detail = 0, clusters = NULL, theta = NULL,
    theta.des = NULL, parfunc = NULL, dparfunc = NULL,
    step = 1, sym = 0, colnames = NULL, dimpar = NULL,
    weights = NULL, score.method = "nlminb",
    same.cens = FALSE, censoring.probs = NULL, silent = 1,
    ...)
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

  \item{parfunc}{R-function to specify any relationship for
  the parameters of the cross-odds-ratio.}

  \item{dparfunc}{R-function to specify the derivative of
  the parfunc function with respect to the parameters.}

  \item{step}{specifies the step size for the
  Newton-Raphson algorithm.}

  \item{sym}{specifies if symmetry is used in the model.}

  \item{colnames}{specifies possible colnames to go with
  the parameters.}

  \item{dimpar}{number of parameters in case where
  functions specify the model.}

  \item{weights}{weights for estimating equations.}

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
functd<-function(theta,t,x) { xt=theta[1]+theta[2]*(t-0.5); return(xt) };
dfunctd<-function(theta,t,x) { return(c(x[1],(t-0.5))) };

#library(compiler) # using the compiler library to speed up things
#f=cmpfun(functd)
#df=cmpfun(dfunctd)
#out4<-cor.cif(add,data=multcif,multcif$cause,
#parfunc=f,dparfunc=df,
#cause1=1,cause2=1,cens.model="KM",clusters=multcif$id,theta=log(2+1),
#theta.des=theta.des,colnames=c("intercept","slope"))
#summary(out4)

###out4<-cor.cif(add,data=multcif,parfunc=functd,dparfunc=dfunctd,
###             cause1=1,cause2=1,cens.model="KM",theta=log(2+1),
###             theta.des=theta.des,colnames=c("intercept","slope"),dimpar=2)
###summary(out4)
## no indication of timetrend in cross-odds-ratio
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

