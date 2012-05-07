\name{twinlm}
\alias{twinlm}
\title{Classic twin model for quantitative traits}
\usage{
  twinlm(formula, data, id, zyg, DZ, DZos, weight = NULL,
    type = c("ace"), twinnum = "twinnum", binary = FALSE,
    keep = weight, estimator = "gaussian", ...)
}
\arguments{
  \item{formula}{Formula specifying effects of covariates
  on the response.}

  \item{data}{\code{data.frame} with one observation pr
  row. In addition a column with the zygosity (DZ or MZ
  given as a factor) of each individual much be specified
  as well as a twin id variable giving a unique pair of
  numbers/factors to each twin pair.}

  \item{id}{The name of the column in the dataset
  containing the twin-id variable.}

  \item{zyg}{The name of the column in the dataset
  containing the zygosity variable.}

  \item{DZ}{Character defining the level in the zyg
  variable corresponding to the dyzogitic twins. If this
  argument is missing, the reference level (i.e. the first
  level) will be interpreted as the dyzogitic twins.}

  \item{DZos}{Optional. Character defining the level in the
  zyg variable corresponding to the oppposite sex dyzogitic
  twins.}

  \item{weight}{Weight matrix if needed by the chosen
  estimator. For use with Inverse Probability Weights}

  \item{type}{Character defining the type of analysis to be
  performed. Should be a subset of "aced" (additive genetic
  factors, common environmental factors, unique
  environmental factors, dominant genetic factors).}

  \item{twinnum}{The name of the column in the dataset
  numbering the twins (1,2). If it does not exist in
  \code{data} it will automatically be created.}

  \item{binary}{If \code{TRUE} a liability model is fitted}

  \item{keep}{Vector of variables from \code{data} that are
  not specified in \code{formula}, to be added to
  data.frame of the SEM}

  \item{estimator}{Choice of estimator/model.}

  \item{...}{Additional arguments parsed on to lower-level
  functions}
}
\value{
  Returns an object of class \code{twinlm}.
}
\description{
  Fits a classical twin model for quantitative traits.
}
\examples{
## Simulate data
set.seed(1)
d <- twinsim(2000,b1=c(1,-1),b2=c(),acde=c(1,1,0,1))
## E(y|z1,z2) = z1 - z2. var(A) = var(C) = var(E) = 1

## E.g to fit the data to an ACE-model without any confounders we simply write
ace <- twinlm(y1 ~ 1, data=d, DZ="DZ", zyg="zyg", id="id")
ace
## An AE-model could be fitted as
ae <- twinlm(y1 ~ 1, data=d, DZ="DZ", zyg="zyg", id="id", type="ae")
## LRT:
compare(ae,ace)
## AIC
AIC(ae)-AIC(ace)
## To adjust for the covariates we simply alter the formula statement
ace2 <- twinlm(y1 ~ x11+x12, data=d, DZ="DZ", zyg="zyg", id="id", type="ace")
## Summary/GOF
summary(ace2)
## An interaction could be analyzed as:
ace3 <- twinlm(y1 ~ x11+x12 + x11:I(x12<0), data=d, DZ="DZ", zyg="zyg", id="id", type="ace")
## Categorical variables are also supported
d2 <- transform(d,x12cat=cut(x12,3,labels=c("Low","Med","High")))
ace4 <- twinlm(y1 ~ x11+x12cat, data=d2, DZ="DZ", zyg="zyg", id="id", type="ace")
## plot the model structure
\donttest{
plot(ace4)
}
}
\author{
  Klaus K. Holst
}
\seealso{
  \code{\link{twinsim}}
}
\keyword{models}
\keyword{regression}

