\name{easy.binomial.twostage}
\alias{easy.binomial.twostage}
\title{Fits two-stage binomial for describing depdendence in binomial data
using marginals that are on logistic form using the binomial.twostage funcion, but
call is different and easier and the data manipulation is build into the function.
Useful in particular for family design data.}
\usage{
  easy.binomial.twostage(margbin = NULL,
    data = sys.parent(), score.method = "nlminb",
    response = "response", id = "id", Nit = 60, detail = 0,
    silent = 1, weights = NULL, control = list(),
    theta = NULL, theta.formula = NULL, desnames = NULL,
    deshelp = 0, var.link = 1, iid = 1, step = 0.5,
    model = "plackett", marginal.p = NULL, strata = NULL,
    max.clust = NULL, se.clusters = NULL)
}
\arguments{
  \item{margbin}{Marginal binomial model}

  \item{data}{data frame}

  \item{response}{name of response variable in data frame}

  \item{id}{name of cluster variable in data frame}

  \item{score.method}{Scoring method}

  \item{Nit}{Number of iterations}

  \item{detail}{Detail for more output for iterations}

  \item{silent}{Debug information}

  \item{weights}{Weights for log-likelihood, can be used
  for each type of outcome in 2x2 tables.}

  \item{control}{Optimization arguments}

  \item{theta}{Starting values for variance components}

  \item{theta.formula}{design for depedence, either formula
  or design function}

  \item{desnames}{names for dependence parameters}

  \item{deshelp}{if 1 then prints out some data sets that
  are used, on on which the design function operates}

  \item{var.link}{Link function for variance}

  \item{iid}{Calculate i.i.d. decomposition}

  \item{step}{Step size}

  \item{model}{model}

  \item{marginal.p}{vector of marginal probabilities}

  \item{strata}{strata for fitting}

  \item{max.clust}{max clusters}

  \item{se.clusters}{clusters for iid decomposition for
  roubst standard errors}
}
\description{
  If clusters contain more than two times, the algoritm
  uses a compososite likelihood based on the pairwise
  bivariate models.
}
\details{
  The reported standard errors are based on the estimated
  information from the likelihood assuming that the
  marginals are known. This gives correct standard errors
  in the case of the plackett distribution (OR model for
  dependence), but incorrect for the clayton-oakes types
  model. The OR model is often known as the ALR model, but
  this function gives correct standard errors and is quite
  a bit quicker.
}
\examples{
data(twinstut)
twinstut0 <- subset(twinstut, tvparnr<2300000)
theta.des <- model.matrix( ~-1+factor(zyg),data=twinstut0)
margbin <- glm(stutter~factor(sex)+age,data=twinstut0,family=binomial())
bin <- binomial.twostage(margbin,data=twinstut0,clusters=twinstut0$tvparnr,
			    theta.des=theta.des,detail=0,
	  score.method="fisher.scoring",Nit=20,step=1.0)
summary(bin)

##bin <- binomial.twostage(stutter~factor(sex)+age,data=twinstut0,clusters=twinstut0$tvparnr,
## 			    theta.des=theta.des,detail=0,
## 	  score.method="fisher.scoring",Nit=20,step=1.0)
## summary(bin)

### easier call with formula call for OR dependence model
twinstut0$binstut <- (twinstut0$stutter=="yes")*1
out <- easy.binomial.twostage(stutter~factor(sex)+age,data=twinstut0,response="binstut",id="tvparnr",
			      theta.formula=~-1+factor(zyg1))
summary(out)

### easier call with flexible design for OR dependence model using functions
desfs <- function(x,num1="zyg1",namesdes=c("mz","dz","os")) c(x[num1]=="dz",x[num1]=="mz",x[num1]=="os")*1
out3 <- easy.binomial.twostage(binstut~factor(sex)+age,
      data=twinstut0,response="binstut",id="tvparnr",
      score.method="fisher.scoring", theta.formula=desfs,desnames=c("mz","dz","os"))
summary(out3)

### Family design with parents and children and dependence parameters for
### mother - father (mf), mother-child and father-child (mb), child-child (bb)
dd <- mets:::sim.bin.fam(1000,beta=0,theta=2)
ddl <- fast.reshape(dd,varying="y",keep="y")

### design function sets up the type of the pair and the related design
desfs <- function(x,num1="num1",num2="num2")
{
	     mf <- (x[num1]=="m")*(x[num2]=="f")*1
         mb <- (x[num1]=="m" | x[num1]=="f")*(x[num2]=="b1" | x[num2]=="b2")*1
         bb <- (x[num1]=="b1")*(x[num2]=="b1" | x[num2]=="b2")*1
         c(mf,mb,bb)
}

out <- easy.binomial.twostage(y~+1,data=ddl,response="y",id="id",
        	       score.method="fisher.scoring",deshelp=0,
		       theta.formula=desfs,desnames=c("mf","mb","bb"))
summary(out)
}
\keyword{binomial}
\keyword{regression}

