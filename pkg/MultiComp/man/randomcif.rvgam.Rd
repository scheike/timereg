\name{Grandom.cif}
\alias{Grandom.cif}
\title{Random effects model for competing risks data} 
\description{
Fits a random effects  model describing the dependence in the cumulative 
incidence curves for subjects within a cluster.  Given the gamma distributed
random effects it is assumed that the cumulative incidence curves 
are indpendent, and that the marginal cumulative incidence curves 
are on additive form
\deqn{
P(T \leq t, cause=1 | x,z) = P_1(t,x,z) = 1- exp( -x^T A(t) - t z^T \beta)
}

We allow a regression structure for the indenpendent gamma distributed 
random effects  and their variances that may depend on cluster covariates.

Design.rv specificies the random effects for each cluster. This is
a matrix of 1's and 0's with dimension n x d.  With d random effects. 
For a cluster with two subjects, we let the design.rv rows be 
 \eqn{v_1} and \eqn{v_2}. 
Such that the random effects for subject 
1 is \deqn{v_1^T (Z_1,...,Z_d)}, for d random effects. Each random effect
has an associated parameter \eqn{(\lambda_1,...,\lambda_d)}. By construction
subjects 1's random effect are Gamma distributed with 
mean \eqn{\lambda_1/v_1^T \lambda}
and variance \eqn{\lambda_1/(v_1^T \lambda)^2}. Note that the random effect 
\eqn{v_1^T (Z_1,...,Z_d)} has mean 1 and variance \eqn{1/(v_1^T \lambda)}.

The parameters \eqn{(\lambda_1,...,\lambda_d)}
are related to the parameters of the model
by a regression construction \eqn{pard} (d x k), that links the \eqn{d} 
\eqn{\lambda} parameters
with the (k) underlying \eqn{\theta} parameters 
\deqn{
\lambda = pard \theta 
}
}
\usage{
Grandom.cif(cif,data=sys.parent(),cause,
causeS=1,cens.code=0,cens.model="KM",Nit=40,detail=0,
clusters=NULL,theta=NULL,theta.des=NULL,step=1,design.rv=NULL,link=0,
notaylor=1,same.cens=FALSE,entry=NULL,trunkp=1)
}
\arguments{
\item{cif}{a model object from the comp.risk function with the 
marginal cumulative incidence of cause2, i.e., the event that is conditioned on, and whose
odds the comparision is made with respect to}
\item{data}{a data.frame with the variables.}
\item{cause}{specifies the causes  related to the death
	times, the value cens.code is the censoring value.}
\item{causeS}{specificies the cause of interst.}
\item{cens.code}{specificies the code for the censoring.}
\item{cens.model}{specified which model to use for the ICPW, KM is Kaplan-Meier alternatively it may be "cox"}
\item{Nit}{number of iterations for Newton-Raphson algorithm.}
\item{detail}{if 0 no details are printed during iterations, if 1 details are given.}
\item{clusters}{specifies the cluster structure.}
\item{theta}{specifies starting values for the cross-odds-ratio parameters of the model.}
\item{design.rv}{specifies a regression design for the random effects.}
\item{link}{specifies a possible exponential link function.}
\item{theta.des}{specifies a regression design for the parameters of the 
random effects, each row most give the deign. Default is independent parameters,
that is identity matrices of dimension equivalent to number of random effects.}
\item{step}{specifies the step size for the Newton-Raphson algorithm.}
\item{notaylor}{if TRUE then the standard errors ignores the uncertainty 
from the marginal cumulative incidence. This will speed things up considerably.} 
\item{same.cens}{if true then censoring within clusters are assumed to be the same variable, default is independent censoring.}
\item{entry}{entryage for delayed entry.}
\item{trunkp}{truncation probabilities, ie P(T > entryage) estimated from other model.}
}
\value{returns an object of type 'cor'. With the following arguments:
\item{theta}{estimate of parameters of model.}
\item{var.theta}{variance for gamma.  }
\item{hess}{the derivative of the used score.}
\item{score}{scores at final stage.}
\item{theta.iid}{matrix of iid decomposition of parametric effects.}
}
\references{
A Semiparametric Random Effects Model for Multivariate Competing Risks Data,
Scheike, Zhang, Sun, Jensen (2010), Biometrika.

Cross odds ratio Modelling of dependence for
Multivariate Competing Risks Data, Scheike and Sun (2009), work in progress.


Scheike, Holst, Hjelmborg (2010), 
}
\author{Thomas Scheike}
\examples{
data(multcif)
multcif$cause[multcif$cause==0] <- 2

addm<-comp.risk(Surv(time,status>0)~const(X)+cluster(id),data=multcif,
               multcif$cause,n.sim=0)

### making group indidcator 
g.des<-data.frame(group2=rep(rbinom(200,1,0.5),rep(2,200)))
theta.des <- model.matrix(~-1+factor(group2),g.des)

out1m<-random.cif(addm,data=multcif,causeS=1,Nit=15,detail=0,
theta=2,theta.des=theta.des,step=1.0)
summary(out1m)

## this model can also be formulated as a random effects model 
## but different parameters
out2m<-Grandom.cif(addm,data=multcif,causeS=1,Nit=10,detail=0,
theta=1/c(4.84,0.79),design.rv=theta.des,step=1.0)
summary(out2m)
1/out2m$theta
out1m$theta

####################################################################
################### AC modelling of twin data ######################
####################################################################
### assume that zygbin gives the zygosity of mono and dizygotic twins
### 0 for mono and 1 for dizygotic twins. We now formulate and AC model
zygbin <- g.des$group2 ## indicator of dizygotic twins

n <- nrow(multcif)
### random effects for each cluster
des.rv <- cbind(theta.des,(zygbin==1)*rep(c(1,0)),(zygbin==1)*rep(c(0,1)),1)
### design making parameters half the variance for dizygotic components
pard <- rbind(c(1,0), c(0.5,0),c(0.5,0), c(0.5,0), c(0,1))
pardes <- matrix(pard,n,10,byrow=TRUE)

outacem <-Grandom.cif(addm,data=multcif,causeS=1,Nit=30,detail=0,
          theta=c(-1.21,2.1),theta.des=pardes,step=1.0,design.rv=des.rv,
          link=0,notaylor=1)
summary(outacem)
### genetic variance is 
outacem$theta[1]/sum(outacem$theta)^2
### best if variances had been positive
}
\keyword{survival}
