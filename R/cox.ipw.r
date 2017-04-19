#' Missing data IPW Cox
#' 
#' Fits an Cox-Aalen survival model with missing data, with glm specification
#' of probability of missingness.
#' 
#' Taylor expansion of Cox's partial likelihood in direction of glm parameters
#' using num-deriv and iid expansion of Cox and glm paramters (lava).
#' 
#' @aliases cox.ipw summary.cox.ipw print.cox.ipw coef.cox.ipw
#' @param survformula a formula object with the response on the left of a '~'
#' operator, and the independent terms on the right as regressors. The response
#' must be a survival object as returned by the `Surv' function.
#' 
#' Adds the prop() wrapper internally for using cox.aalen function for fitting
#' Cox model.
#' @param glmformula formula for "being" observed, that is not missing.
#' @param d data frame.
#' @param max.clust number of clusters in iid approximation. Default is all.
#' @param ipw.se if TRUE computes standard errors based on iid decompositon of
#' cox and glm model, thus should be asymptotically correct.
#' @param tie.seed if there are ties these are broken, and to get same break
#' the seed must be the same. Recommend to break them prior to entering the
#' program.
#' @return returns an object of type "cox.aalen". With the following arguments:
#' \item{iid}{iid decomposition.} \item{coef}{missing data estiamtes for
#' weighted cox. } \item{var}{robust pointwise variances estimates.  }
#' \item{se}{robust pointwise variances estimates.  } \item{se.naive}{estimate
#' of parametric components of model.  } \item{ties}{list of ties and times
#' with random noise to break ties.} \item{cox}{output from weighted cox
#' model.}
#' @author Thomas Scheike
#' @references Paik et al.
#' @keywords survival
#' @examples
#' 
#' 
#' ### fit <- cox.ipw(Surv(time,status)~X+Z,obs~Z+X+time+status,data=d,ipw.se=TRUE)
#' ### summary(fit)
#' 
#' 
##' @export
cox.ipw <- function(survformula,glmformula,d=sys.parent(),max.clust=NULL,ipw.se=FALSE,tie.seed=100)
{ ## {{{ 
  ggl <- glm(glmformula,family='binomial',data=d)
  mat <-  model.matrix(glmformula,data=d);
  glmcovs <- attr(ggl$terms,"term.labels")
  d$ppp <- predict(ggl,type='response')

###  d1 <- d[,survcovs]
###  dcc <- na.omit(d)
  ## {{{ checking and breaking ties
  ties <- FALSE
  survtimes <- all.vars(update(survformula,.~0))
  if (length(survtimes)==2) {itime <- 1; time2 <- d[,survtimes[1]]; status <- d[,survtimes[2]]; }
  if (length(survtimes)==3) {itime <- 2;  time2 <- d[,survtimes[2]]; status <- d[,survtimes[3]]; }
  jtimes <- time2[status==1]; 
  dupli <- duplicated(jtimes)
  if (sum(dupli)>0) {
	  set.seed(tie.seed)
	  jtimes[dupli] <- jtimes[dupli]+runif(sum(dupli))*0.01
	  time2[status==1] <- jtimes
	  d[,survtimes[itime]] <- time2
	  ties <- TRUE
  }
  ## }}}

  dcc <- d[ggl$y==1,]
  ppp <- dcc$ppp
  timeregsurvformula <- timereg.formula(survformula)
  udca <- cox.aalen(timeregsurvformula,data=dcc,weights=1/ppp,n.sim=0,max.clust=max.clust)  

  ### iid of beta for Cox model 
  coxiid <- udca$gamma.iid

if (ipw.se==TRUE)  { ## {{{ 
###requireNamespace("lava"); requireNamespace("NumDeriv"); 
glmiid <-   lava::iid(ggl)
mat <- mat[ggl$y==1,]
par <- coef(ggl)

coxalpha <- function(par)
{ ## {{{ 
  rr <- mat %*% par
  pw <- c(exp(rr)/(1+exp(rr)))
  assign("pw",pw,envir=environment(survformula))
###  if (coxph==FALSE) 
  ud <- cox.aalen(timeregsurvformula,data=dcc,weights=1/pw,beta=udca$gamma,Nit=1,n.sim=0,robust=0)  
###  else { ud <- coxph(survformula,data=dcc,weights=1/pw,iter.max=1,init=udca$gamma)  
###	 ud <- coxph.detail(ud,data=dcc)
###  }
  ud$score
} ## }}} 

DU <-  numDeriv::jacobian(coxalpha,par,)
IDU <-  udca$D2linv %*% DU 
alphaiid <-t( IDU %*% t(glmiid))
###
iidfull <- alphaiid
###
iidfull[ggl$y==1,] <- coxiid + alphaiid[ggl$y==1,]
###
var2 <- t(iidfull) %*% iidfull
se <- cbind(diag(var2)^.5); colnames(se) <- "se"
} else { iidfull <- NULL; var2 <- NULL; se <- NULL} ## }}} 

se.naive=coef(udca)[,3,drop=FALSE]; colnames(se.naive) <- "se.naive"

res <- list(iid=iidfull,coef=udca$gamma,var=var2,se=se,se.naive=se.naive,ties=list(ties=ties,time2=time2,cox=udca))
class(res) <- "cox.ipw"
return(res)
} ## }}} 

##' @export
summary.cox.ipw <- function(object,digits=3,...)
{
	tval <- object$coef/object$se
	pval <- 2*(1-pnorm(abs(tval)))
       	res <- cbind(object$coef,object$se,object$se.naive,pval)
	colnames(res) <- c("coef","se","se.naive","pval")

return(res)
}

##' @export
coef.cox.ipw<- function(object,digits=3,...)
{
summary.cox.ipw(object)
}

##' @export
print.cox.ipw  <-  function(x,...)
{ 
summary.cox.ipw(x)
}

