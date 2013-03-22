\name{easy.binomial.twostage}
\alias{easy.binomial.twostage}
\title{Fits two-stage binomial for describing depdendence in binomial data}
\usage{
  easy.binomial.twostage(margbin = NULL,
    data = sys.parent(), score.method = "nlminb",
    response = "response", id = "id", type = NULL,
    Nit = 60, detail = 0, silent = 1, weights = NULL,
    control = list(), theta = NULL, theta.formula = NULL,
    desnames = NULL, deshelp = 0, var.link = 1, iid = 1,
    step = 0.5, model = "plackett", marginal.p = NULL,
    strata = NULL, max.clust = NULL, se.clusters = NULL)
}
\description{
  Fits two-stage binomial for describing depdendence in
  binomial data
}
\examples{
data(twinstut)
theta.des <- model.matrix( ~-1+factor(zyg),data=twinstut)
margbin <- glm(stutter~factor(sex)+age,data=twinstut,family=binomial())
bin <- binomial.twostage(margbin,data=twinstut,clusters=twinstut$tvparnr,
			    theta.des=theta.des,detail=0,
	  score.method="fisher.scoring",Nit=20,step=1.0)
summary(bin)

bin <- binomial.twostage(stutter~factor(sex)+age,data=twinstut,clusters=twinstut$tvparnr,
			    theta.des=theta.des,detail=0,
	  score.method="fisher.scoring",Nit=20,step=1.0)
summary(bin)

twinstut$binstut <- (twinstut$stutter=="yes")*1
out <- easy.binomial.twostage(stutter~factor(sex)+age,data=twinstut,response="binstut",id="tvparnr",
			      theta.formula=~-1+factor(zyg1))
summary(out)

desfs <- function(x,num1="zyg1",num2="num2",namesdes=c("mz","dz","os")) c(x[num1]=="mz",x[num1]=="dz",x[num1]=="os")*1
out3 <- easy.binomial.twostage(binstut~factor(sex)+age,
      data=twinstut,response="binstut",id="tvparnr",type="zyg",
      score.method="fisher.scoring", theta.formula=desfs,desnames=c("mz","dz","os"))
summary(out3)
}

