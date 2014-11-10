cox.ipw <- function(survformula,glmformula,d=sys.parent(),max.clust=NULL,ipw.se=FALSE)
{ ## {{{ 
  ggl <- glm(glmformula,family='binomial',data=d)
  glmcovs <- attr(ggl$terms,"term.labels")
  d$ppp <- predict(ggl,type='response')

###  d1 <- d[,survcovs]
###  dcc <- na.omit(d)
  dcc <- d[ggl$y==1,]
  ppp <- dcc$ppp
  udca <- cox.aalen(survformula,data=dcc,weights=1/ppp,n.sim=0,max.clust=max.clust)  

  ### iid of beta for Cox model 
###  coxiid <- t(solve(udca$D2linv) %*% t(udca$gamma.iid))
  coxiid <- udca$gamma.iid

if (ipw.se==TRUE)  { ## {{{ 
if (!require(numDeriv)) stop("numDeriv needed\n")
if (!require(lava))     stop("lava needed\n")
glmiid <-   lava::iid(ggl)
mat <-  model.matrix(glmformula,data=dcc);
par <- coef(ggl)

coxalpha <- function(par)
{ ## {{{ 
  rr <- mat %*% par
  pw <- c(exp(rr)/(1+exp(rr)))
  assign("pw",pw,envir=environment(survformula))
  ud <- cox.aalen(survformula,data=dcc,weights=1/pw,beta=udca$gamma,Nit=1,n.sim=0,robust=0)  
  ud$score
} ## }}} 

DU <-  numDeriv::jacobian(coxalpha,par)
IDU <-  udca$D2linv %*% DU 
alphaiid <-t( IDU %*% t(glmiid))
###
iidfull <- alphaiid
###
iidfull[ggl$y==1,] <- coxiid - alphaiid[ggl$y==1,]
iidfull <-t( udca$D2linv %*% t(iidfull))
###
var2 <- t(iidfull) %*% iidfull
se <- cbind(diag(var2)^.5); colnames(se) <- "se"
} else { var2 <- NULL; se <- NULL} ## }}} 

se.naive=coef(udca)[,3,drop=FALSE]; colnames(se.naive) <- "se.naive"

res <- list(iid=iidfull,coef=udca$gamma,var=var2,se=se,se.naive=se.naive)
return(res)
} ## }}} 

