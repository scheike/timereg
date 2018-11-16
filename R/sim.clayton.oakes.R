##' Simulate observations from the Clayton-Oakes copula model with
##' piecewise constant marginals.
##'
##' @title Simulate from the Clayton-Oakes frailty model
##' @param K Number of clusters
##' @param n Number of observations in each cluster
##' @param eta variance
##' @param beta Effect (log hazard ratio) of covariate
##' @param stoptime Stopping time
##' @param lam constant hazard 
##' @param left Left truncation
##' @param pairleft pairwise (1) left truncation or individual (0)
##' @param trunc.prob Truncation probability 
##' @param same if 1 then left-truncation is same also for univariate truncation
##' @author Thomas Scheike and Klaus K. Holst
##' @aliases simClaytonOakes
##' @export
simClaytonOakes <- function(K,n,eta,beta,stoptime,lam=1,left=0,pairleft=0,trunc.prob=0.5,same=0)  ## {{{ 
{
  ## K antal clustre, n=antal i clustre
  ###	K <- 100; n=2; stoptime=2; eta=1/2; beta=0; lam=0.5;left=0.5; trunc.prob=0.5; pairleft=0; same=0
  ### change such that eta is variance 
###  eta <- 1/eta
  x<-array(c(runif(n*K),rep(0,n*K),rbinom(n*K,1,0.5)),dim=c(K,n,3))
  C<-matrix(stoptime,K,n);
  Gam1 <-matrix(rgamma(K,eta)/eta,K,n)
  temp<-eta*log(-log(1-x[,,1])/(eta*Gam1*lam)+1)*exp(-beta*x[,,3])
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
       if (same==0) lefttime <- rexp(n*K)*left
       if (same==1) lefttime <- rep(rexp(K)*left,each=n)
###       lefttime <- runif(K)*(stoptime-left)
       if (same==0) left <- rbinom(n*K,1,trunc.prob) ## not trunation times!
       if (same==1) left <- rep(rbinom(K,1,trunc.prob),each=n)
       lefttime <- lefttime*left
       trunk <- ud[,1] > lefttime
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
  names(ud)<-c("time","status","x","cluster","mintime","lefttime","truncated")
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
simClaytonOakesWei <- function(K,n,eta,beta,stoptime,weiscale=1,weishape=2,left=0,pairleft=0)
{ ## {{{ 
###cat(" not quite \n"); 
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

  colnames(ud)<-c("time","status","x","cluster","mintime","lefttime","truncated")
  ud <- data.frame(ud)
  return(ud)
} ## }}} 

##' @export
simClaytonOakes.twin.ace <- function(K,varg,varc,beta,stoptime,Cvar=0,left=0,pairleft=0,trunc.prob=0.5,lam0=1)  ## {{{ 
{
  ## K antal clustre, n=antal i clustre
  n <- 2 # twins with ace structure
  #change parametrization 
  sumpar <- sum(varg+varc)
  varg <- varg/sumpar^2; 
  varc <- varc/sumpar^2
  x<-array(c(runif(n*K),rep(0,n*K),rbinom(n*K,1,0.5)),dim=c(K,n,3))
  if (Cvar==0) C<-matrix(stoptime,K,n) else C<-matrix(Cvar*runif(K*n)*stoptime,K,n) 
  ### total variance of gene and env. 
  ###  random effects with 
  ###  means varg/(varg+varc) and variances varg/(varg+varc)^2
  etao <- eta <- varc+varg
  if (etao==0) eta <- 1
  Gams1 <-cbind( rgamma(K,varg)/eta, rgamma(K,varg*0.5)/eta, rgamma(K,varg*0.5)/eta, rgamma(K,varg*0.5)/eta, rgamma(K,varc)/eta )
###  print(apply(Gams1,2,mean)); print(apply(Gams1,2,var))
  mz <- c(rep(1,K/2),rep(0,K/2)); dz <- 1-mz;
  mzrv <-  Gams1[,1]+Gams1[,5]           ### shared gene + env 
  dzrv1 <- Gams1[,2]+Gams1[,3]+Gams1[,5] ### 0.5 shared gene + 0.5 non-shared + env 
  dzrv2 <- Gams1[,2]+Gams1[,4]+Gams1[,5] ### 0.5 shared gene + 0.5 non-shared + env 
  Gam1 <- cbind(mz*mzrv+dz*dzrv1,mz*mzrv+dz*dzrv2)
###  print(apply(Gam1,2,mean)); print(apply(Gam1,2,var))
  Gam1[Gam1==0] <- 1 ## to work also under independence 
###  print(mean(mzrv)); print(mean(dzrv1)); print(mean(dzrv2)); 
###  print(var(mzrv));  print(var(dzrv1));  print(var(dzrv2)); 
  temp<-eta*log(-log(1-x[,,1])/(eta*Gam1)+1)*exp(-beta*x[,,3])/lam0
  if (etao==0) temp <- matrix(rexp(n*K),K,n)*exp(-beta*x[,,3])/lam0
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

names(ud)<-c("time","status","x","cluster","zyg","mintime","lefttime","truncated")
return(ud)
} ## }}} 

##' @export
simClaytonOakes.family.ace <- function(K,varg,varc,beta,stoptime,lam0=0.5,Cvar=0,left=0,pairleft=0,trunc.prob=0.5)  ## {{{ 
{
  ## K antal clustre (families), n=antal i clustre
  n <- 4 # twins with ace structure
  x<- array(c(runif(n*K),rep(0,n*K),rbinom(n*K,1,0.5)),dim=c(K,n,3))
  if (Cvar==0) C<-matrix(stoptime,K,n) else C<-matrix(Cvar*runif(K*n)*stoptime,K,n) 
  sumpar <- sum(varg+varc)
  varg <- varg/sumpar^2; 
  varc <- varc/sumpar^2
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
  Gam1 <- cbind(mother,father,child1,child2)
###  print(apply(Gam1,2,mean)); print(apply(Gam1,2,var))
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
  if (pairleft==1) ud <- cbind(ud,type,rep(minstime,each=n),rep(lefttime,each=n),rep(trunk,each=n))
  else ud <- cbind(ud,type,rep(minstime,each=n),lefttime,trunk)

names(ud)<-c("time","status","x","cluster","type","mintime","lefttime","truncated")
return(ud)
} ## }}} 


##' @export
simCompete.twin.ace <- function(K,varg,varc,beta,stoptime,lam0=c(0.2,0.3),
		Cvar=0,left=0,pairleft=0,trunc.prob=0.5,overall=1,all.sum=1)  ## {{{ 
{
  ## K antal clustre, n=antal i clustre
  n=2 # twins with ace structure
  sumpar <- sum(varg+varc)
  varg <- varg/sumpar^2; 
  varc <- varc/sumpar^2
  ## length(lam0) competing risk with constant hazards lam0
  x<-array(c(runif(n*K),rep(0,n*K),rbinom(n*K,1,0.5)),dim=c(K,n,3))
  if (Cvar==0) C<-matrix(stoptime,K,n) else C<-matrix(Cvar*runif(K*n)*stoptime,K,n) 
  ### total variance of gene and env. 
  ### one for each cause and one shared (across causes)
  ###  random effects with 
  ###  means varg/(varg+varc) and variances varg/(varg+varc)^2
  if (length(varc)==1) varc  <- rep(varc,length(lam0)+overall)
  if (length(varg)==1) varg  <- rep(varg,length(lam0)+overall)
  eta <- varc+varg
  etat <- sum(eta)
  ### total variance for each cause + overall
  nc <- length(lam0); 
###  etat <- sum(eta[1:nc])
###  print(etat)
###  print(varc)
###  print(varg)
###  print("MZ shared variance"); print(eta[1]/eta[1]^2); 
###  print("DZ shared variance"); 
###  print(c(mdz,vdz,vdz/mdz^2)); 

  mz <- c(rep(1,K/2),rep(0,K/2)); dz <- 1-mz;
  ### ace overall 
  if (overall==1) {
    varcl <- varc[nc+1]; vargl <- varg[nc+1]
    if (all.sum==1) etal <-  etat else etal <- vargl+varcl
    Gams1 <-cbind(
       rgamma(K,vargl)/etal, rgamma(K,vargl*0.5)/etal, rgamma(K,vargl*0.5)/etal, rgamma(K,vargl*0.5)/etal,
       rgamma(K,varcl)/etal )
###  ex1 <- Gams1[,6]
###  ex2 <- Gams1[,7]
  mzrv <-  Gams1[,1]+          Gams1[,5] ### shared gene + env 
  dzrv1 <- Gams1[,2]+Gams1[,3]+Gams1[,5] 
  dzrv2 <- Gams1[,2]+Gams1[,4]+Gams1[,5]
  Gamoa <- cbind(mz*mzrv+dz*dzrv1,mz*mzrv+dz*dzrv2)
  } else Gamoa <- 0

  temp1 <- matrix(0,K,length(lam0))
  temp2 <- matrix(0,K,length(lam0))
  Gamm <- c()
  for (i in 1:nc)
  {
  varcl <- varc[i]; vargl <- varg[i]
  if (all.sum==1) etal <-  etat else etal <- vargl+varcl
  Gams1 <-cbind(
      rgamma(K,vargl)/etal, 
      rgamma(K,vargl*0.5)/etal, rgamma(K,vargl*0.5)/etal, rgamma(K,vargl*0.5)/etal,
      rgamma(K,varcl)/etal )
###  
  mzrv <-  Gams1[,1]+Gams1[,5]     ### shared gene + env 
  dzrv1 <- Gams1[,2]+Gams1[,3]+Gams1[,5] 
  dzrv2 <- Gams1[,2]+Gams1[,4]+Gams1[,5]
  Gam1 <- cbind(mz*mzrv+dz*dzrv1,mz*mzrv+dz*dzrv2)
  Gam1 <- Gam1+Gamoa
  Gamm <- cbind(Gamm,Gam1)
### ## {{{ 
###  mean(mzrv)
###  var(mzrv)
###  var(mzrv[mz==1])
###  mean(dzrv1)
###  mean(dzrv2)
###  var(dzrv1)
###  var(dzrv2)
###  apply(Gam1,2,mean); apply(Gam1,2,var)
###  apply(Gam1[mz==1,],2,mean); apply(Gam1[mz==0,],2,mean)
###  var(Gam1[mz==1,]); var(Gam1[mz==0,])
###  shdz <- Gams1[,2]+Gams1[,5]
###  mean(shdz)
###  ###1.25/etat
###  var(shdz)
###  mean(shdz/0.83)
###  var(shdz/0.83)
###  ## }}}
  ttemp<-matrix(rexp(2*K),K,2)/(Gam1*exp(beta*x[,,3])*lam0[i])
  temp1[,i] <- ttemp[,1]
  temp2[,i] <- ttemp[,2]
  }
###  print(cov(Gamm))
###  temp0 <- cbind( temp1, temp2)
###  print(cov(temp0))
  temp <- cbind( apply(temp1,1,min), apply(temp2,1,min))
  cause1 <- apply(temp1,1,which.min)
  cause2 <- apply(temp2,1,which.min)
###  for (zyg in c(0,1)) 
###  for (i in 1:2) for (j in 1:2) {
###	  med <- (i==cause1) & (j==cause2) & mz==zyg
###	  datl <- temp[med,]
###	  dato <- temp0[med,]
###	  print(c(zyg,i,j))
###	  print(cor(datl))
###	  print(c(zyg,i,j))
###	  print(cor(dato))
###  }
###
  x[,,2]<- ifelse(temp<=C,1,0)*cbind(cause1,cause2);
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

names(ud)<-c("time","status","x","cluster","zyg","mintime","lefttime","truncated")
return(ud)
} ## }}} 

##' @export
simCompete.simple <- function(K,varr,beta,stoptime,lam0=c(0.2,0.3),
	Cvar=0,left=0,pairleft=0,trunc.prob=0.5,overall=1,all.sum=1)  ## {{{ 
{
  ## K antal clustre, n=antal i clustre
  n=2 # twins with ace structure
  ## length(lam0) competing risk with constant hazards lam0
  x<-array(c(runif(n*K),rep(0,n*K),rbinom(n*K,1,0.5)),dim=c(K,n,3))
  if (Cvar==0) C<-matrix(stoptime,K,n) else C<-matrix(Cvar*runif(K*n)*stoptime,K,n) 
  ## variance and mean of additve gamma via paramenters 
  sp <- sum(varr)
  partheta <- varr/sp^2 

  eta <- partheta
  etat <- sum(eta)
  ### total variance for each cause + overall
  nc <- length(lam0); 
  print(eta)
  print(etat)

  mz <- c(rep(1,K/2),rep(0,K/2)); dz <- 1-mz;
  ### ace overall 
  if (overall==1) {
  ###    if (all.sum==1) etal <-  etat else etal <- vargl+varcl
    etal <- etat
    Gams1 <-cbind(rgamma(K,eta[nc+1])/etal)
  Gamoa <- Gams1
  } else Gamoa <- 0

###  print(apply(Gamoa,2,mean))
###  print(apply(Gamoa,2,var))

  temp1 <- matrix(0,K,length(lam0))
  temp2 <- matrix(0,K,length(lam0))
  for (i in 1:nc)
  {
	  etal <- etat
	  Gams1 <-cbind(rgamma(K,eta[i])/etal)
	  Gam1 <- Gams1+Gamoa
	  Gam1 <- cbind(Gam1,Gam1)
###  print("_________________")
###  print(apply(Gam1,2,mean))
###  print(apply(Gam1,2,var))
	  occ <- (1:nc)[-i]
	  for (j in  occ) {
		  print(eta[j])
	     Gamo <- cbind(rgamma(K,eta[j])/etal,rgamma(K,eta[j])/etal)
###             print(apply(Gamo,2,mean))
###             print(apply(Gamo,2,var))
	     Gam1 <- Gam1+Gamo
	  }
  print(apply(Gam1,2,mean))
  print(apply(Gam1,2,var))
	  ttemp<-matrix(rexp(2*K),K,2)/(Gam1*exp(beta*x[,,3])*lam0[i])
	  temp1[,i] <- ttemp[,1]
	  temp2[,i] <- ttemp[,2]
  }
  temp <- cbind( apply(temp1,1,min), apply(temp2,1,min))
  cause1 <- apply(temp1,1,which.min)
  cause2 <- apply(temp2,1,which.min)
###
  x[,,2]<- ifelse(temp<=C,1,0)*cbind(cause1,cause2);
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

names(ud)<-c("time","status","x","cluster","zyg","mintime","lefttime","truncated")
return(ud)
} ## }}} 

##' @export
simFrailty.simple <- function(K,varr,beta,stoptime,lam0=c(0.2),
		Cvar=0,left=0,pairleft=0,trunc.prob=0.5,overall=1,all.sum=NULL)  ## {{{ 
{
  n=2 
  ## length(lam0) competing risk with constant hazards lam0
  x<-array(c(runif(n*K),rep(0,n*K),rbinom(n*K,1,0.5)),dim=c(K,n,3))
  if (Cvar==0) C<-matrix(stoptime,K,n) else C<-matrix(Cvar*runif(K*n)*stoptime,K,n) 
  if (length(varr)==1) varr  <- rep(varr,length(lam0)+overall)
  eta <- varr
  etat <- sum(varr)
  if (!is.null(all.sum)) etat <- all.sum 
  varr <- varr/etat 
  ### total variance for each cause + overall
  nc <- length(lam0); 

  mz <- c(rep(1,K/2),rep(0,K/2)); dz <- 1-mz;
  if (overall==1) {
    etal <- etat
    Gams1 <-cbind(rgamma(K,varr[nc+1])/etal)
    Gamoa <- Gams1
  } else Gamoa <- 0

  print(etat); print(varr)
  temp1 <- matrix(0,K,length(lam0))
  temp2 <- matrix(0,K,length(lam0))
  for (i in 1:nc)
  {
  etal <- etat
  Gams1 <-cbind(rgamma(K,varr[i])/etal)
  Gam1 <- Gams1+Gamoa
  Gam1 <- cbind(Gam1,Gam1)
  print(i); 
  print(apply(Gam1,2,mean)); print(apply(Gam1,2,var)); 
  print(apply(Gams1,2,mean)); print(apply(Gams1,2,var)); 
  print(apply(as.matrix(Gamoa),2,mean)); print(apply(as.matrix(Gamoa),2,var)); 
  ttemp<-matrix(rexp(2*K),K,2)/(Gam1*exp(beta*x[,,3])*lam0[i])
  temp1[,i] <- ttemp[,1]
  temp2[,i] <- ttemp[,2]
  }
  temp <- cbind( apply(temp1,1,min), apply(temp2,1,min))
  cause1 <- apply(temp1,1,which.min)
  cause2 <- apply(temp2,1,which.min)
###
  x[,,2]<- ifelse(temp<=C,1,0)*cbind(cause1,cause2);
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

names(ud)<-c("time","status","x","cluster","zyg","mintime","lefttime","truncated")
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
  sumpar <- sum(parg+parc)
  parg <- parg/sumpar^2; 
  parc <- parc/sumpar^2
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

  ## Silence false R CMD CHECK warnings:
  V11 <- V12 <- V21 <- V22 <- NULL
    
  gams.pair <- fast.reshape(Gam1,id="id")
  gams.pair <-  transform(gams.pair,
              kendall = ((V11-V12)*(V21-V22))/((V11+V12)*(V21+V22))
  )

  kendall <- gams.pair$kendall
  mz.kendall <- mean(kendall[gams.pair$mz1==1])
  dz.kendall <- mean(kendall[gams.pair$mz1==0])

  return(list(mz.kendall=mz.kendall,dz.kendall=dz.kendall))
} ## }}} 

##' @export
kendall.normal.twin.ace <- function(parg,parc,K=10000)  ## {{{ 
{
  ## K antal clustre, n=antal i clustre
  ### total variance of gene and env. 
  ###  K <- 10; varg <- 1; varc <- 1; 
  ###  random effects with 
  ###  means varg/(varg+varc) and variances varg/(varg+varc)^2
  K <- K*2
  Gams1 <-cbind( parg^.5*rnorm(K,parg), (parg*0.5)^.5*rnorm(K), 
                 (parg*0.5)^.5*rnorm(K), (parg*0.5)^.5*rnorm(K), parc^.5*rnorm(K) )
  mz <- c(rep(1,K/2),rep(0,K/2)); 
  dz <- 1-mz;
  id <- rep(1:(K/2),each=2)
  mzrv <-  Gams1[,1]+Gams1[,5]           ### shared gene + env 
  dzrv1 <- Gams1[,2]+Gams1[,3]+Gams1[,5] ### 0.5 shared gene + 0.5 non-shared + env 
  dzrv2 <- Gams1[,2]+Gams1[,4]+Gams1[,5] ### 0.5 shared gene + 0.5 non-shared + env 
  Gam1 <- cbind(mz*mzrv+dz*dzrv1,mz*mzrv+dz*dzrv2)
  Gam1 <- data.frame(cbind(exp(Gam1),mz,id))

  ## Silence false R CMD CHECK warnings:
  V11 <- V12 <- V21 <- V22 <- NULL

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

