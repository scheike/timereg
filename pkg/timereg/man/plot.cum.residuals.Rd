\name{plot.cum.residuals}
\alias{plot.cum.residuals}
\title{Plots cumulative residuals}
\description{This function plots the output from the 
cumulative residuals function "cum.residuals". 
The cumulative residuals are  compared with the performance of similar
processes under the model.
}
\usage{
\method{plot}{cum.residuals}(x,pointwise.ci=1,hw.ci=0,sim.ci=0,
robust=1, specific.comps=FALSE,level=0.05,start.time=0,stop.time=0,
add.to.plot=FALSE,mains=TRUE,xlab="Time",
ylab ="Cumulative Residuals",ylim=NULL,score=0,conf.band=FALSE,...)
}
\arguments{
\item{x}{the output from the "cum.residuals" function.}
\item{pointwise.ci}{if >1 pointwise confidence intervals are
plotted with lty=pointwise.ci
}
\item{hw.ci}{if >1 Hall-Wellner confidence bands are
plotted with lty=hw.ci. Only 95\% bands can be constructed.
}
\item{sim.ci}{if >1 simulation based confidence bands are
plotted with lty=sim.ci. These confidence bands are robust to
non-martingale behaviour.
}
\item{robust}{if "1" robust standard errors are used to estimate standard 
error of estimate, otherwise martingale based estimate are used.
}
\item{specific.comps}{all components of the model is plotted by default,
but a list of components may be specified, for example first and third
"c(1,3)".
}
\item{level}{gives the significance level. Default is 0.05.
}
\item{start.time}{
start of observation period where estimates are plotted. Default is 0.}
\item{stop.time}{end of period where estimates are plotted.
Estimates thus plotted from [start.time, max.time].}
\item{add.to.plot}{to add to an already existing plot. Default is "FALSE".
}
\item{mains}{add names of covariates as titles to plots.}
\item{xlab}{label for x-axis. Default is "Time".}
\item{ylab}{label for y-axis. Default is "Cumulative Residuals".}
\item{ylim}{limits for y-axis.}
\item{score}{if '0' plots related to modelmatrix are specified, thus
  resulting in grouped residuals, if '1' plots for modelmatrix but with
  random realizations under model, if '2' plots residuals versus
  continuous covariates of model with random realizations under the
  model.}
\item{conf.band}{makes simulation based confidence bands for the test processes 
under the 0 based on variance of these processes limits for y-axis. These will
give additional information of whether the observed cumulative residuals are 
extreme or not when based on a variance weighted test.}
\item{...}{unused arguments - for S3 compatibility}

}
\references{
Martinussen and Scheike, Dynamic Regression Models for Survival Data,
Springer (2006).
}
\author{Thomas Scheike}
\examples{
# see cum.residuals for examples 
}
\keyword{survival}
