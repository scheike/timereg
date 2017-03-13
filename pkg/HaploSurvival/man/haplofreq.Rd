\name{haplo.freqs}
\alias{haplo.freqs}
\alias{print.haplo.freqs}
\title{Fit Hardy-Weinberg equilibrium model for haplotype frequencies based on genotype data}
\description{
Fit MLE of Hardy-Weinberg equilibrium model for haplotype frequencies based on genotype data 

Input is genotype data. 

The haplotype model assumes Hardy-Weinberg equilibrium 
\deqn{ P(H=(h_k,h_j)) = \pi_k pi_j }
with K different haplotypes, where we parametrize the frequencies as
\deqn{ 
\pi_k= \frac{ \exp(\theta_k) }{\exp(\theta_1)+...+\exp(\theta_{K-1})+1}
}
We allow a regression design on the haplotype parameters to reduce the 
dimensionality 
\deqn{ \theta = G \alpha }
where G is matrix of dimension (K-1) x (d-1). The \eqn{\alpha} and \eqn{\theta}
version both have the same baseline that is a single haplotype.  
}
\usage{
haplo.freqs(geno.type, geno.setup=NULL, Nit = 10,detail = 0,
haplo.freq = NULL, step = 1,lev.marq=1,min.lev.marq=0, 
haplo.design=NULL,haplo.baseline=NULL,alpha=NULL)
}
\arguments{
\item{geno.type}{a data.frame with the variables.}
\item{geno.setup}{setup form genotype data.}
\item{Nit}{number of iterations for Newton-Raphson algorithm.}
\item{detail}{if 0 no details is printed during iterations, if 1 details are given.}
\item{haplo.freq}{haplo type frequencies for starting values or for fixed values.}
\item{step}{step size for iteration steps.}
\item{lev.marq}{starting value for Levenberg-Marquart algorithm.}
\item{min.lev.marq}{stopping value for Levenberg-Marquart algorithm after initial steps, min.lev.marq=0 is
stopping value for Levenberg-Marquart algorithm after initial steps, min.lev.marq=0 is
Newton-Raphson steps, that are needed to get correct standard errors.}
\item{haplo.design}{design for haplo parameters, default is diagonal of size K-1}
\item{haplo.baseline}{gives the name of the geno.setup that should be used as baseline, default
is last haplotype.}
\item{alpha}{starting value of parameters for  haplo.design}
}
\value{
returns an object of type "haplo.freq". With the following arguments:
\item{haplo.alpha}{MLE of alpha parameters for haplo.design.}
\item{haplo.pars}{MLE of haplotype parameters.}
\item{haplo.freq}{MLE of haplotype frequencies.}
\item{var.haplo.alpha}{variance matrix for alpha}
\item{D2linv}{second derivative of log-likelihood.}
\item{score}{score for likelihood.}
\item{haplo.design}{design that relates haplotype parameters to alpha parameters.}
\item{alpha.iid}{iid decomposition of alpha parameters, in this case score contributions 
for each subject for score of for log-likelihood.}
\item{freq.names}{names of haplotypes}
}
\references{
Scheike, Martinussen and Silver,  Haplotype effects for survival data, Submitted.

}
\author{Thomas Scheike and Jeremy Silver}
\examples{
n<-200
g<-matrix(rbinom(200*4,1,0.5),n,4); # simple genotypes

gs<-geno.setup(g); 

hgs<-haplo.freqs(g,geno.setup=gs)
summary(hgs); 
hgs$haplo.freq         # haplo-type frequencies
hgs$haplo.pars         #  related parameters 
hgs$haplo.alpha        #  related parameters  haplo.pars = X alpha
hgs$var.haplo.alpha    #  variance of alpha parameters 

# structured haplo-frequency model, baseline is set to "0,0"

gs<-geno.setup(g,haplo.baseline="0,0"); 

X<-matrix(1,3,1); 
hgs<-haplo.freqs(g,geno.setup=gs,haplo.design=X)
summary(hgs); 
hgs$haplo.freq         # haplo-type frequencies
hgs$haplo.pars         # related parameters 
hgs$haplo.alpha        # related parameters  haplo.pars = X alpha
hgs$var.haplo.alpha    # variance of alpha parameters 

}
\keyword{survival}
