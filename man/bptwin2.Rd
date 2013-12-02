\name{bptwin2}
\alias{bptwin2}
\title{Liability model for twin data}
\usage{
  bptwin2(formula, data, id, zyg, DZ, group = NULL,
    num = NULL, weight = NULL,
    biweight = function(x) 1/min(x), strata = NULL,
    messages = 1, control = list(trace = 0), type = "ace",
    eqmean = TRUE, pairsonly = FALSE, samecens = TRUE,
    allmarg = samecens & !is.null(weight), stderr = TRUE,
    robustvar = TRUE, p, indiv = FALSE, constrain,
    bound = FALSE, varlink, ...)
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
  variable corresponding to the dyzogitic twins.}

  \item{group}{Optional. Variable name defining group for
  interaction analysis (e.g., gender)}

  \item{num}{Optional twin number variable}

  \item{weight}{Weight matrix if needed by the chosen
  estimator. For use with Inverse Probability Weights}

  \item{biweight}{Function defining the bivariate weight in
  each cluster}

  \item{strata}{Strata}

  \item{messages}{Control amount of messages shown}

  \item{control}{Control argument parsed on to the
  optimization routine. Starting values may be parsed as
  '\code{start}'.}

  \item{type}{Character defining the type of analysis to be
  performed. Should be a subset of "aced" (additive genetic
  factors, common environmental factors, unique
  environmental factors, dominant genetic factors).}

  \item{eqmean}{Equal means (with type="u")?}

  \item{pairsonly}{Include complete pairs only?}

  \item{stderr}{Should standard errors be calculated?}

  \item{robustvar}{If TRUE robust (sandwich) variance
  estimates of the variance are used}

  \item{p}{Parameter vector p in which to evaluate
  log-Likelihood and score function}

  \item{indiv}{If TRUE the score and log-Likelihood
  contribution of each twin-pair}

  \item{constrain}{Development argument}

  \item{samecens}{Same censoring}

  \item{allmarg}{Should all marginal terms be included}

  \item{bound}{Development argument}

  \item{varlink}{Link function for variance parameters}

  \item{...}{Additional arguments to lower level functions}
}
\description{
  Liability-threshold model for twin data
}
\examples{
\donttest{
data(twinstut)
b0 <- bptwin(stutter~sex,
             data=droplevels(subset(twinstut,zyg\%in\%c("mz","dz"))),
             id="tvparnr",zyg="zyg",DZ="dz",type="ae")
summary(b0)
}
}
\author{
  Klaus K. Holst
}

