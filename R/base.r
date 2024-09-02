coefBase<- function(object, digits=3, d2logl=0,ci=1,alpha=0.05) { ## {{{
  res <- cbind(object$gamma,
               diag(object$var.gamma)^0.5,
               diag(object$robvar.gamma)^0.5)
  if (d2logl==1) res<-cbind(res,diag(object$D2linv)^.5)
  wald <- object$gamma/diag(object$robvar.gamma)^0.5
  waldp <- (1 - pnorm(abs(wald))) * 2
  res <- signif(as.matrix(cbind(res, wald, waldp)),digits=digits)
  if (d2logl==1) colnames(res) <- c("Coef.", "SE", "Robust SE","D2log(L)^-1","z","P-val") else colnames(res) <- c("Coef.", "SE", "Robust SE", "z", "P-val")
  if (ci==1) {
     slower <- paste("lower",signif(100*alpha/2,2),"%",sep="")
     supper <- paste("upper",signif(100*(1-alpha/2),3),"%",sep="")
     res <- signif(cbind(res,res[,1]+qnorm(alpha/2)*res[,2],res[,1]+qnorm(1-alpha/2)*res[,2]),digits=digits); 
     nn <- ncol(res); colnames(res)[(nn-1):nn] <- c(slower,supper)
  }
###  prmatrix(signif(res, digits))
  return(res)
} ## }}}

coefcox <- function(object, digits=3, d2logl=0,ci=1,alpha=0.05) { ## {{{
  res <- cbind(object$coef, diag(object$var)^0.5)
  wald <- object$coef/diag(object$var)^0.5
  waldp <- (1 - pnorm(abs(wald))) * 2
  res <- signif(as.matrix(cbind(res, wald, waldp)),digits=digits)
  colnames(res) <- c("Coef.", "SE", "z", "P-val")
  if (ci==1) {
     slower <- paste("lower",signif(100*alpha/2,2),"%",sep="")
     supper <- paste("upper",signif(100*(1-alpha/2),3),"%",sep="")
     res <- signif(cbind(res,res[,1]+qnorm(alpha/2)*res[,2],res[,1]+qnorm(1-alpha/2)*res[,2]),digits=digits); 
     nn <- ncol(res); colnames(res)[(nn-1):nn] <- c(slower,supper)
  }
###  prmatrix(signif(res, digits))
  return(res)
} ## }}}


#' Makes wald test
#' 
#' Makes wald test, either by contrast matrix or testing components to 0. Can
#' also specify the regression coefficients and the variance matrix.  Also
#' makes confidence intervals of the defined contrasts.  Reads coefficientes
#' and variances from timereg and coxph objects.
#' 
#' 
#' @param object timereg object
#' @param coef estimates from some model
#' @param Sigma variance of estimates
#' @param vcov same as Sigma but more standard in other functions
#' @param contrast contrast matrix for testing
#' @param coef.null which indeces to test to 0
#' @param null mean of null, 0 by default
#' @param print.coef print the coefficients of the linear combinations.
#' @param alpha significance level for CI for linear combinations of
#' coefficients.
#' @author Thomas Scheike
#' @keywords survival
#' @examples
#' 
#' data(sTRACE)
#' # Fits Cox model 
#' out<-cox.aalen(Surv(time,status==9)~prop(age)+prop(sex)+
#' prop(vf)+prop(chf)+prop(diabetes),data=sTRACE,n.sim=0)
#' 
#' wald.test(out,coef.null=c(1,2,3))
#' ### test age=sex   vf=chf
#' wald.test(out,contrast=rbind(c(1,-1,0,0,0),c(0,0,1,-1,0)))
#' 
#' ### now same with direct specifation of estimates and variance
#' wald.test(coef=out$gamma,Sigma=out$var.gamma,coef.null=c(1,2,3))
#' wald.test(coef=out$gamma,Sigma=out$robvar.gamma,coef.null=c(1,2,3))
#' ### test age=sex   vf=chf
#' wald.test(coef=out$gamma,Sigma=out$var.gamma,
#' 	  contrast=rbind(c(1,-1,0,0,0),c(0,0,1,-1,0)))
#' 
#' @export
wald.test <- function(object=NULL,coef=NULL,Sigma=NULL,vcov=NULL,contrast,coef.null=NULL,null=NULL,print.coef=TRUE,alpha=0.05)
{ ## {{{

  coefs <- NULL
  if (inherits(object,"coxph"))  {coef <-  matrix(coef(object),ncol=1); Sigma=object$var;}
  if (inherits(object,"phreg"))  {coef <-  matrix(c(coef(object)),ncol=1); Sigma=vcov(object);}
  if (inherits(object,"cox.aalen"))  {coef <- object$gamma; Sigma=object$var.gamma;}
  if (is.null(Sigma)) {
     if (inherits(object,c("cor","twostage"))) Sigma <- object$var.theta else Sigma <- object$var.gamma;
  }
  if (!is.null(object)) {
     if (inherits(object,c("cor","twostage"))) coefs <- object$theta else coefs <- object$gamma;
  } 
  if (is.null(coefs)) coefs <- coef(object)
  if (!is.null(coef)) coefs <- coef ## else stop("No estimates given \n"); 
  nl <- length(coefs)

  if (missing(contrast)) {
      contrast <- rep(1,length(coefs))
      contrast <- diag(1,nl);
###      namesc <- names(c(coefs))
  } 
  if (!is.null(coef.null)) {
	  contrast <- c()
	  for (i in coef.null) 
      contrast <- rbind(contrast,c((1:nl)==i)*1)
  }
  if (missing(null)) null <- 0

  if (is.null(Sigma)) Sigma <- vcov
  if (is.null(Sigma)) Sigma <- vcov(object)

 ### Wald test
 B <- contrast
 p <- coefs
 if (is.vector(B)) { B <- rbind(B); colnames(B) <- names(contrast) }

 varBp <- B%*%Sigma%*%t(B)
 seBp <- diag(varBp)^.5
 lin.comb <- B %*% p
 Q <- t(B%*%p-null)%*%solve(varBp)%*%(B%*%p-null)
 zvals <- lin.comb/seBp
 pvals <- 2*(1-pnorm(abs(zvals)))
 coef.out <- cbind(lin.comb,seBp,lin.comb+seBp*qnorm(alpha/2),lin.comb+seBp*qnorm(1-alpha/2),pvals)
 colnames(coef.out) <- c("lin.comb","se","lower","upper","pval")
 if (print.coef) prmatrix(coef.out)

 df <- qr(B)$rank; names(df) <- "df"
 attributes(Q) <- NULL; names(Q) <- "chisq";
 pQ <- ifelse(df==0,NA,1-pchisq(Q,df))
 method = "Wald test";
 ##    hypothesis <-
 res <- list(##data.name=hypothesis,
  statistic = Q, parameter = df, p.value=pQ, method = method,coef.out=coef.out,varBp=varBp,lin.comb=lin.comb)
  class(res) <- "htest"
  attributes(res)$B <- B
return(res)
} ## }}}


timetest<-function(object,digits=3,hyp.label="p-value H_0:constant effect",out=0)
{  ## {{{
  cat("Test for nonparametric terms \n")
  if (is.null(object$conf.band)==TRUE)  mtest<-FALSE else mtest<-TRUE;
  if (mtest==FALSE) cat("Test not computed, sim=0 \n\n")
  if (mtest==TRUE) {
  test0<-cbind(object$obs.testBeq0,object$pval.testBeq0)
  testC<-cbind(object$obs.testBeqC,object$pval.testBeqC)
  colnames(test0)<-c("Supremum-test of significance","p-value H_0: B(t)=0")
  colnames(testC)<-c("      Kolmogorov-Smirnov test",hyp.label)
  if (is.null(object$obs.testBeqC.is)!=TRUE)  {
  testCis<-cbind(object$obs.testBeqC.is,object$pval.testBeqC.is)
  colnames(testCis) <-
                   c("        Cramer von Mises test",hyp.label)
  }
  cat("\n")
  cat("Test for non-significant effects \n")
  prmatrix(signif(test0,digits))
  cat("\n")
  cat("Test for time invariant effects \n")
  prmatrix(signif(testC,digits))
  if (is.null(object$obs.testBeqC.is)!=TRUE)  prmatrix(signif(testCis,digits))
  cat("\n")
  if (out==1) return(cbind(test0,testC)); 
}
} ## }}}

is.diag <- function(m)
{ ## {{{
p <- nrow(m)
adiag <- min(diag(m)*1)
if (adiag==0) ud <- FALSE else ud <- TRUE
dm <- diag(p); diag(dm) <- diag(m); 
ndiag <- sum(abs(c(m - dm)))
if (ndiag>0.0000001) ud <- FALSE;
return(ud)
} ## }}}


risk.index <- function(start,stop,id,times)
{ ## {{{
n <- length(start)
nstop <- length(stop)
if (n!=nstop) stop("start and stop not of same length\n"); 
if (is.null(id)) id <- 1:n
nid <- length(id)
if (n!=nid) stop("id and start not of same length\n"); 

nt <- length(times)

 nclust <- .C("atriskindex",
	as.double(start), as.double(stop), as.integer(id), as.integer(n),
	as.double(times), as.integer(nt), as.integer(rep(0,nt)),as.integer(rep(0,nt*n)),PACKAGE="timereg")

  nrisk <- nclust[[7]]
  riskindex <- matrix(nclust[[8]],nt,n)

out <- list(nrisk=nrisk,riskindex=riskindex)
} ## }}}

namematrix<-function(mat,names)
{ colnames(mat)<-names; rownames(mat)<-names; return(mat); }
nameestimate<-function(mat,names)
{ colnames(mat)<-"estimate"; rownames(mat)<-names; return(mat); }

aalen.des<-function(formula=formula(data),data=parent.frame(),model="aalen")
{ ## {{{
  call <- match.call(); 
  m <- match.call(expand.dots=FALSE); 
  m$model<-NULL
  special <- c("cluster","prop","const")
  Terms <- if(missing(data)) terms(formula,special) else terms(formula, special, data=data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  mt <- attr(m, "terms")
  intercept<-attr(mt, "intercept")
  Y <- model.extract(m,"response")

  des<-read.design(m,Terms,model=model)
  X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px;
  pz<-des$pz;
  covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ;
  clusters<-des$clusters;

  if (attr(m[, 1], "type") == "right") {
    type<-"right"; 
    status <- m[, 1][, "status"];
    time2  <- m[, 1][, "time"]; time   <- rep(0,length(time2)); 
  } else if (attr(m[, 1], "type") == "counting") {
    type<-"counting"; 
    time   <- m[, 1][,1]; time2  <- m[, 1][,2]; status <- m[, 1][,3];
  } else { stop("only right-censored or counting processes data") } 
return(list(type=type,time=time,time2=time2,status=status,
 X=X,Z=Z,px=px,pz=pz,npar=npar,
 covnamesX=covnamesX,covnamesZ=covnamesZ,clusters=clusters))
} ## }}}

#' @export
names2formula <- function(formula,names)
{ ## {{{ 
covs <- paste("~",names[1],sep="")
if (length(names)>=2) 
for (i in 2:length(names))
covs <- c(paste(covs,paste("+",names[i],sep="")))
covs <- as.formula(covs)
formula <- update(formula,covs)

return(formula)
} ## }}} 


#' @export
timereg.formula <- function(formula,propterms=NULL,special="prop")
{ ## {{{ 
vars <- all.vars(update(formula,0~.))
nonpropvars <- NULL
if (!is.null(propterms)) { nonpropvars <- vars[-propterms]; 
                           vars <- vars[propterms] 
}

covs <- paste("~",special,"(",vars[1],")",sep="")
if (length(vars)>=2) 
for (i in 2:length(vars))
covs <- c(paste(covs,paste("+",special,"(",vars[i],")",sep="")))

if (!is.null(nonpropvars))
{ 
for (i in 1:length(nonpropvars))
covs <- c(paste(covs,paste("+",nonpropvars[i],sep="")))
}
covs <- as.formula(covs)

formula <- update(formula,covs)
return(formula)
} ## }}} 

