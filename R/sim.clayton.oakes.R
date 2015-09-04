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
simClaytonOakes.twin.ace <- function(K,varg,varc,beta,stoptime,Cvar=0,left=0,pairleft=0,trunc.prob=0.5)  ## {{{ 
{
  ## K antal clustre, n=antal i clustre
  n=2 # twins with ace structure
  x<-array(c(runif(n*K),rep(0,n*K),rbinom(n*K,1,0.5)),dim=c(K,n,3))
  if (Cvar==0) C<-matrix(stoptime,K,n) else C<-matrix(cvar*runif(K*n)*stoptime,K,n) 
  ### total variance of gene and env. 
  ###  random effects with 
  ###  means varg/(varg+varc) and variances varg/(varg+varc)^2
  eta <- varc+varg
  Gams1 <-cbind(
       rgamma(K,varg)/eta,
       rgamma(K,varg*0.5)/eta, rgamma(K,varg*0.5)/eta, rgamma(K,varg*0.5)/eta,
       rgamma(K,varc)/eta )
  mz <- c(rep(1,K/2),rep(0,K/2)); dz <- 1-mz;
  mzrv <-  Gams1[,1]+Gams1[,5]           ### shared gene + env 
  dzrv1 <- Gams1[,2]+Gams1[,3]+Gams1[,5] ### 0.5 shared gene + 0.5 non-shared + env 
  dzrv2 <- Gams1[,2]+Gams1[,4]+Gams1[,5] ### 0.5 shared gene + 0.5 non-shared + env 
  Gam1 <- cbind(mz*mzrv+dz*dzrv1,mz*mzrv+dz*dzrv2)
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

##' @export
simClaytonOakes.family.ace <- function(K,varg,varc,beta,stoptime,lam0=0.5,Cvar=0,left=0,pairleft=0,trunc.prob=0.5)  ## {{{ 
{
  ## K antal clustre (families), n=antal i clustre
  n=4 # twins with ace structure
  x<-array(c(runif(n*K),rep(0,n*K),rbinom(n*K,1,0.5)),dim=c(K,n,3))
  if (Cvar==0) C<-matrix(stoptime,K,n) else C<-matrix(cvar*runif(K*n)*stoptime,K,n) 
  ### total variance of gene and env. 
  ###  random effects with 
  ###  means varg/(varg+varc) and variances varg/(varg+varc)^2
  eta <- varc+varg
  ### mother and father share environment
  ### children share half the genes with mother and father and environment 
  mother.g <-  cbind(rgamma(K,varg*0.25)/eta, rgamma(K,varg*0.25)/eta, rgamma(K,varg*0.25)/eta, rgamma(K,varg*0.25)/eta)
  father.g <-  cbind(rgamma(K,varg*0.25)/eta, rgamma(K,varg*0.25)/eta, rgamma(K,varg*0.25)/eta, rgamma(K,varg*0.25)/eta)
  env <- rgamma(K,varc)/eta 
  mother <- apply(mother.g,1,sum)+env
  father <- apply(father.g,1,sum)+env
  child1 <- apply(cbind(mother.g[,c(1,2)],father.g[,c(1,2)]),1,sum) + env
  child2 <- apply(cbind(mother.g[,c(1,3)],father.g[,c(1,3)]),1,sum) + env
###  Gams1 <-cbind(
###       rgamma(K,varg*0.25)/eta, rgamma(K,varg*0.25)/eta, rgamma(K,varg*0.25)/eta, rgamma(K,varg*0.25)/eta,
###       rgamma(K,varg*0.5)/eta, rgamma(K,varg*0.5)/eta, rgamma(K,varg*0.5)/eta,
###       rgamma(K,varc)/eta )
###  mz <- c(rep(1,K/2),rep(0,K/2)); dz <- 1-mz;
###  mzrv <-  Gams1[,1]+Gams1[,5]           ### shared gene + env 
###  dzrv1 <- Gams1[,2]+Gams1[,3]+Gams1[,5] ### 0.5 shared gene + 0.5 non-shared + env 
###  dzrv2 <- Gams1[,2]+Gams1[,4]+Gams1[,5] ### 0.5 shared gene + 0.5 non-shared + env 
###  Gam1 <- cbind(mz*mzrv+dz*dzrv1,mz*mzrv+dz*dzrv2)
  Gam1 <- cbind(mother,father,child1,child2)
  temp<-eta*log(-log(1-x[,,1])/(eta*Gam1)+1)*exp(-beta*x[,,3])/lam0
  x[,,2]<-ifelse(temp<=C,1,0);
  x[,,1]<-pmin(temp,C)
  minstime <- apply(x[,,1],1,min)  
  ud <- as.data.frame(cbind(apply(x,3,t),rep(1:K,each=n)))  
  type <- rep(c("mother","father","child","child"),K)

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
  else ud <- cbind(ud,type,rep(minstime,each=n),lefttime,trunk)

names(ud)<-c("time","status","x1","cluster","type","mintime","lefttime","truncated")
return(ud)
} ## }}} 


##' @export
kendall.ClaytonOakes.twin.ace <- function(parg,parc,K=10000)  ## {{{ 
{
  ## K antal clustre, n=antal i clustre
  ### total variance of gene and env. 
  ###  K <- 10; varg <- 1; varc <- 1; 
  ###  random effects with 
  ###  means varg/(varg+varc) and variances varg/(varg+varc)^2
  K <- K*2
  eta <- parc+parg
  Gams1 <-cbind( rgamma(K,parg)/eta, rgamma(K,parg*0.5)/eta, 
                 rgamma(K,parg*0.5)/eta, rgamma(K,parg*0.5)/eta, 
		 rgamma(K,parc)/eta )
  mz <- c(rep(1,K/2),rep(0,K/2)); 
  dz <- 1-mz;
  id <- rep(1:(K/2),each=2)
  mzrv <-  Gams1[,1]+Gams1[,5]           ### shared gene + env 
  dzrv1 <- Gams1[,2]+Gams1[,3]+Gams1[,5] ### 0.5 shared gene + 0.5 non-shared + env 
  dzrv2 <- Gams1[,2]+Gams1[,4]+Gams1[,5] ### 0.5 shared gene + 0.5 non-shared + env 
  Gam1 <- cbind(mz*mzrv+dz*dzrv1,mz*mzrv+dz*dzrv2)
  Gam1 <- data.frame(cbind(Gam1,mz,id))

  gams.pair <- fast.reshape(Gam1,id="id")
  gams.pair <-  transform(gams.pair,
              kendall = ((V11-V12)*(V21-V22))/((V11+V12)*(V21+V22))
  )

  kendall <- gams.pair$kendall
  mz <- gams.pair$mz1
  mz.kendall <- mean(kendall[mz==1])
  dz.kendall <- mean(kendall[mz==0])

  return(list(mz.kendall=mz.kendall,dz.kendall=dz.kendall))
} ## }}} 
##
###kendall.ClaytonOakes.twin.ace(2,0)

## sim.clayton <- function(n=100,K=2,eta=0.5,beta,...) {
##     m <- lvm(T~x)
##     rates <- c(0.3,0.5); cuts <- c(0,5)
##     distribution(m,~) <- coxExponential.lvm(rate=rates,timecut=cuts)    
## }

