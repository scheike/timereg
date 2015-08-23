##' Simulate observations from the Clayton-Oakes copula model with
##' piecewise constant marginals.
##'
##' @title Simulate from the Clayton-Oakes frailty model
##' @param K Number of clusters
##' @param n Number of observations in each cluster
##' @param eta 1/variance
##' @param beta Effect (log hazard ratio) of covariate
##' @param stoptime Stopping time
##' @param left Left truncation
##' @param pairleft pairwise (1) left truncation or individual (0)
##' @param trunc.prob Truncation probability 
##' @author Klaus K. Holst
##' @export
simClaytonOakes <- function(K,n,eta,beta,stoptime,left=0,pairleft=0,trunc.prob=0.5)  ## {{{ 
{
  ## K antal clustre, n=antal i clustre
  x<-array(c(runif(n*K),rep(0,n*K),rbinom(n*K,1,0.5)),dim=c(K,n,3))
  C<-matrix(stoptime,K,n);
  Gam1 <-matrix(rgamma(K,eta)/eta,K,n)
  temp<-eta*log(-log(1-x[,,1])/(eta*Gam1)+1)*exp(-beta*x[,,3])
  x[,,2]<-ifelse(temp<=C,1,0);
  x[,,1]<-pmin(temp,C)
  minstime <- apply(x[,,1],1,min)  
  ud <- as.data.frame(cbind(apply(x,3,t),rep(1:K,each=n)))  

if (left>0) {
     if (pairleft==1) {
     lefttime <- runif(K)*(stoptime-left)
     left <- rbinom(K,1,trunc.prob) ## not trunation times!
     lefttime <- apply(cbind(lefttime*left,3),1,min)
       trunk <- (lefttime > minstime)
       medleft <- rep(trunk,each=n)
     } else {
###       lefttime <- rexp(n*K)*left
       lefttime <- runif(K)*(stoptime-left)
       left <- rbinom(n*K,1,trunc.prob) ## not trunation times!
       lefttime <- apply(cbind(lefttime*left,3),1,min)
       trunk <- (lefttime > ud[,1])
       medleft <- trunk
     }
  } else { lefttime <- trunk <- rep(0,K);}
  if (pairleft==1) ud <- cbind(ud,rep(minstime,each=n),rep(lefttime,each=n),rep(trunk,each=n))
  else ud <- cbind(ud,rep(minstime,each=n),lefttime,trunk)

###  if (left>0) {
###    lefttime <- rexp(K)*left
###    left <- rbinom(K,1,0.5) ## not trunation times!
###    lefttime <- apply(cbind(lefttime*left,3),1,min)
###    trunk <- (lefttime > minstime)
###    medleft <- rep(trunk,each=n)
###  } else { lefttime <- trunk <- rep(0,K);}
###
###  ud <- cbind(ud,rep(minstime,each=n),rep(lefttime,each=n),rep(trunk,each=n))
  names(ud)<-c("time","status","x1","cluster","mintime","lefttime","truncated")
  return(ud)
} ## }}} 

##' Simulate observations from the Clayton-Oakes copula model with
##' Weibull type baseline and Cox marginals.
##'
##' @title Simulate from the Clayton-Oakes frailty model
##' @param K Number of clusters
##' @param n Number of observations in each cluster
##' @param eta 1/variance
##' @param beta Effect (log hazard ratio) of covariate
##' @param stoptime Stopping time
##' @param weiscale weibull scale parameter 
##' @param weishape weibull shape parameter 
##' @param left Left truncation
##' @param pairleft pairwise (1) left truncation or individual (0)
##' @author Klaus K. Holst 
##' @export
simClaytonOakesWei <- function(K,n,eta,beta,stoptime,
	       weiscale=1,weishape=2,left=0,pairleft=0)
{ ## {{{ 
cat(" not quite \n"); 
## K antal clustre, n=antal i clustre
### K=10; n=2; eta=1; beta=0.3; stoptime=3; lam=0.5; 
### weigamma=2; left=0; pairleft=0
 X <- rbinom(n*K,1,0.5)
 C<-rep(stoptime,n*K);
 Gam1 <-rep(rgamma(K,eta),each=n)
 temp <- rexp(K*n)
### temp <- rweibull(n*K,weishape,scale=weiscale)/(exp(X*beta)*Gam1)
 temp<- (eta*log(eta*temp/(eta*Gam1)+1)/(exp(beta*X)*weiscale^weishape))^{1/weishape}
 status<- ifelse(temp<=C,1,0);
 temp <-   pmin(temp,C)
 xt <- matrix(temp,n,K)
 minstime <- apply(xt,2,min)  
 id=rep(1:K,each=n)
  ud <- cbind(temp,status,X,id)
if (left>0) {
     if (pairleft==1) {
     lefttime <- runif(K)*(stoptime-left)
     left <- rbinom(K,1,0.5) ## not trunation times!
     lefttime <- apply(cbind(lefttime*left,3),1,min)
       trunk <- (lefttime > minstime)
       medleft <- rep(trunk,each=n)
     } else {
###       lefttime <- rexp(n*K)*left
       lefttime <- runif(K)*(stoptime-left)
       left <- rbinom(n*K,1,0.5) ## not trunation times!
       lefttime <- apply(cbind(lefttime*left,3),1,min)
       trunk <- (lefttime > ud[,1])
       medleft <- trunk
     }
  } else { lefttime <- trunk <- rep(0,K);}
  if (pairleft==1) ud <- cbind(ud,rep(minstime,each=n),rep(lefttime,each=n),rep(trunk,each=n))
  else ud <- cbind(ud,rep(minstime,each=n),lefttime,trunk)

  colnames(ud)<-c("time","status","x1","cluster","mintime","lefttime","truncated")
  ud <- data.frame(ud)
  return(ud)
} ## }}} 

##' @export
simClaytonOakes.twin.ace <- function(K,varg,varc,beta,stoptime,pmz=0.5,
				     left=0,pairleft=0,trunc.prob=0.5)  ## {{{ 
{
  ## K antal clustre, n=antal i clustre
  n=2 # twins with ace structure
  x<-array(c(runif(n*K),rep(0,n*K),rbinom(n*K,1,0.5)),dim=c(K,n,3))
  C<-matrix(stoptime,K,n);
  eta <- varc+varg
###  Gam1 <-matrix(rgamma(K,eta)/eta,K,n)
  Gams1 <-cbind(
       rgamma(K,varg)/eta,
       rgamma(K,varg*0.5)/eta, rgamma(K,varg*0.5)/eta, rgamma(K,varg*0.5)/eta,
       rgamma(K,varc)/eta )
###  mz <- rbinom(K,1,pmz); dz <- 1-mz;
  mz <- c(rep(1,K/2),rep(0,K/2)); dz <- 1-mz;
  mzrv <-  Gams1[,1]+Gams1[,5]           ### shared gene + env 
  dzrv1 <- Gams1[,2]+Gams1[,3]+Gams1[,5] ### 0.5 shared gene + 0.5 non-shared + env 
  dzrv2 <- Gams1[,2]+Gams1[,4]+Gams1[,5] ### 0.5 shared gene + 0.5 non-shared + env 
###  mean(mzrv)
###  var(mzrv)
###  mean(dzrv1)
###  var(dzrv1)
###
  Gam1 <- cbind(mz*mzrv+dz*dzrv1,mz*mzrv+dz*dzrv2)
###  apply(Gam1,2,mean)
###  apply(Gam1,2,var)
###mean(gams)
###var(gams)
###cor(Gam1[mz==1,])
###cor(Gam1[dz==1,])
  temp<-eta*log(-log(1-x[,,1])/(eta*Gam1)+1)*exp(-beta*x[,,3])
  x[,,2]<-ifelse(temp<=C,1,0);
  x[,,1]<-pmin(temp,C)
  minstime <- apply(x[,,1],1,min)  
  ud <- as.data.frame(cbind(apply(x,3,t),rep(1:K,each=n)))  
  zyg <- c(rep("MZ",K),rep("DZ",K))

if (left>0) { ## {{{ 
     if (pairleft==1) {
     lefttime <- runif(K)*(stoptime-left)
     left <- rbinom(K,1,trunc.prob) ## not trunation times!
     lefttime <- apply(cbind(lefttime*left,3),1,min)
       trunk <- (lefttime > minstime)
       medleft <- rep(trunk,each=n)
     } else {
###       lefttime <- rexp(n*K)*left
       lefttime <- runif(K)*(stoptime-left)
       left <- rbinom(n*K,1,trunc.prob) ## not trunation times!
       lefttime <- apply(cbind(lefttime*left,3),1,min)
       trunk <- (lefttime > ud[,1])
       medleft <- trunk
     }
  } else { lefttime <- trunk <- rep(0,K);} ## }}} 
  if (pairleft==1) ud <- cbind(ud,zyg,rep(minstime,each=n),rep(lefttime,each=n),rep(trunk,each=n))
  else ud <- cbind(ud,zyg,rep(minstime,each=n),lefttime,trunk)

names(ud)<-c("time","status","x1","cluster","zyg","mintime","lefttime","truncated")
return(ud)
} ## }}} 

###library(mets)
###res <- c()
###for (i in 1:100)
###{
###print(i)
###data <- simClaytonOakes.twin.ace(5000,1,1,0,3)
###names(data)
###table(data$zyg)
######
###out <- polygen.design(data,id="cluster")
###pardes <- out$pardes
###des.rv <- out$des.rv
###tail(des.rv)
######
###aa <- aalen(Surv(time,status)~+1,data=data,robust=0)
###ts <- twostage(aa,data=data,clusters=data$cluster,detail=0,
###	       theta=c(2,1),var.link=0,step=1,
###	       random.design=des.rv,theta.des=pardes)
###summary(ts)
###res <- rbind(res,c(ts$theta,diag(ts$var.theta)^.5))
###}
###apply(res,2,mean)
###apply(res,2,sd)

## sim.clayton <- function(n=100,K=2,eta=0.5,beta,...) {
##     m <- lvm(T~x)
##     rates <- c(0.3,0.5); cuts <- c(0,5)
##     distribution(m,~) <- coxExponential.lvm(rate=rates,timecut=cuts)    
## }

