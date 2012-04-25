##' .. content for description (no empty lines) ..
##'
##' .. content for details ..
##' @title Simulate from the Clayton-Oakes frailty model
##' @param K 
##' @param n 
##' @param eta 
##' @param beta 
##' @param stoptime 
##' @param left 
##' @author Klaus K. Holst
##' @export
simClaytonOakes <- function(K,n,eta,beta,stoptime,left=0) {
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
    lefttime <- rexp(K)*left
    left <- rbinom(K,1,0.5) ## not trunation times!
    lefttime <- apply(cbind(lefttime*left,3),1,min)
    trunk <- (lefttime > minstime)
    medleft <- rep(trunk,each=n)
  } else { lefttime <- trunk <- rep(0,K);}

  ud <- cbind(ud,rep(minstime,each=n),rep(lefttime,each=n),rep(trunk,each=n))
  names(ud)<-c("time","status","x1","cluster","mintime","lefttime","truncated")
  return(ud)
}
