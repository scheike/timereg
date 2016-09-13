\name{cox.ipw}
\alias{cox.ipw}
\alias{summary.cox.ipw}
\alias{print.cox.ipw}
\alias{coef.cox.ipw}
\title{Missing data IPW Cox }
\description{
Fits an Cox-Aalen survival model with missing data, with glm specification of probability of missingness. 
}
\usage{
cox.ipw(survformula,glmformula,d=sys.parent(),max.clust=NULL,ipw.se=FALSE,tie.seed=100)
}
\arguments{
\item{survformula}{
  a formula object with the response on the left of a '~' operator, and
  the independent terms on the right as regressors. The response must be
  a survival object as returned by the `Surv' function.  
  
  Adds the prop() wrapper internally for using cox.aalen function for fitting Cox model.
}
\item{glmformula}{ formula for "being" observed, that is not missing. }
\item{d}{data frame.}
\item{max.clust}{number of clusters in iid approximation. Default is all. }
\item{ipw.se}{if TRUE computes standard errors based on iid decompositon of cox and glm model, thus should be 
asymptotically correct.}
\item{tie.seed}{if there are ties these are broken, and to get same break the seed must be the same. Recommend to 
break them prior to entering the program. }
}
\details{
	Taylor expansion of Cox's partial likelihood in direction of glm parameters using num-deriv
	and iid expansion of Cox and glm paramters (lava). 
}
\value{
returns an object of type "cox.aalen". With the following arguments:
\item{iid}{iid decomposition.}
\item{coef}{missing data estiamtes for weighted cox. }
\item{var}{robust pointwise variances estimates.  }
\item{se}{robust pointwise variances estimates.  }
\item{se.naive}{estimate of parametric components of model.  }
\item{ties}{list of ties and times with random noise to break ties.}
\item{cox}{output from weighted cox model.}
}
\references{
Paik et al. 
}
\author{Thomas Scheike}
\examples{

### fit <- cox.ipw(Surv(time,status)~X+Z,obs~Z+X+time+status,data=d,ipw.se=TRUE)
### summary(fit)

}
\keyword{survival}
