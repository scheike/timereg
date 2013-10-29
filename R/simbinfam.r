
##' @export 
sim.bin.fam <- function(n,beta=0.0,rhopp=0.1,rhomb=0.7,rhofb=0.1,rhobb=0.7) { ## {{{ 
xc <- runif(n)*0.5
xm <- rbinom(n,1,0.5+xc); 
xf <- rbinom(n,1,0.5+xc); 
xb1 <- rbinom(n,1,0.3+xc); 
xb2 <- rbinom(n,1,0.3+xc); 
###
rn <- matrix(rnorm(n*4),n,4)
corm <- matrix( c(1,rhopp,rhomb,rhomb, rhopp,1,rhofb,rhofb, rhomb,rhofb,1,rhobb, rhomb,rhofb,rhobb,1),4,4)
rnn <- t( corm %*% t(rn))
zm <- exp(rnn[,1]); zf <- exp(rnn[,2]); zb1 <- exp(rnn[,3]); zb2 <- exp(rnn[,4]); 
pm <- exp(0.5+xm*beta+zm)
pf <- exp(0.5+xf*beta+zf)
pf <- pf/(1+pf)
pm <- pm/(1+pm)
pb1 <- exp(0.5+xb1*beta+zb1)
pb1 <- pb1/(1+pb1)
pb2 <- exp(0.5+xb2*beta+zb2)
pb2 <- pb2/(1+pb2)
ym <- rbinom(n,1,pm)
yf <- rbinom(n,1,pf)
yb1 <- rbinom(n,1,pb1)
yb2 <- rbinom(n,1,pb2)
#
agem <- 20+runif(n)*10
ageb1 <- 5+runif(n)*10
data.frame(agem=agem,agef=agem+3+rnorm(n)*2,
	   ageb1=ageb1,ageb2=ageb1+1+runif(n)*3,xm=xm,xf=xf,xb1=xb1,xb2=xb2,ym=ym,yf=yf,yb1=yb1,yb2=yb2,id=1:n)
} ## }}} 

