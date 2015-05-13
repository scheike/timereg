\name{plot.aalen}
\alias{plot.aalen}
\alias{plot.cox.aalen}
\alias{plot.timecox}
\alias{plot.prop.excess}
\title{Plots estimates and test-processes}
\description{This function plots the non-parametric cumulative
estimates for the additive risk model or the test-processes for 
the hypothesis of time-varying effects with 
re-sampled processes under the null.
}
\usage{
\method{plot}{aalen}(x,pointwise.ci=1,hw.ci=0,sim.ci=0,robust=0,
specific.comps=FALSE,level=0.05, start.time=0,stop.time=0,add.to.plot=FALSE,
mains=TRUE,xlab="Time",ylab="Cumulative coefficients",score=FALSE,...)
}
\arguments{
\item{x}{the output from the "aalen" function.}
\item{pointwise.ci}{if >1 pointwise confidence intervals are
plotted with lty=pointwise.ci
}
\item{hw.ci}{if >1 Hall-Wellner confidence bands are
plotted with lty=hw.ci. Only 0.95 \% bands can be constructed.
}
\item{sim.ci}{if >1 simulation based confidence bands are
plotted with lty=sim.ci. These confidence bands are robust to
non-martingale behaviour.
}
\item{robust}{robust standard errors are used to estimate standard 
error of estimate, otherwise martingale based standard errors are
used.
}
\item{specific.comps}{all components of the model is plotted by default,
but a list of components may be specified, for example first and third
"c(1,3)".}
\item{level}{gives the significance level.}
\item{start.time}{start of observation period where estimates are plotted.}
\item{stop.time}{end of period where estimates are plotted.
Estimates thus plotted from [start.time, max.time].}
\item{add.to.plot}{to add to an already existing plot.}
\item{mains}{add names of covariates as titles to plots.}
\item{xlab}{label for x-axis.}
\item{ylab}{label for y-axis.}
\item{score}{to plot test processes for test of time-varying effects
  along with 50 random realization under the null-hypothesis.}
\item{...}{unused arguments - for S3 compatibility}
}
\references{
Martinussen and Scheike, Dynamic Regression models for Survival Data,
Springer (2006). 
}
\author{Thomas Scheike}
\examples{
# see help(aalen) 
}
\keyword{survival}
