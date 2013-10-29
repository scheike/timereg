sim.bin.fam <- function(n,beta=0.0,theta=1,lam1=1,lam2=1) { ## {{{ 
xm <- rbinom(n,1,0.5); xf <- rbinom(n,1,0.5); 
xb1 <- rbinom(n,1,0.5); xb2 <- rbinom(n,1,0.5); 
###
zf <- rgamma(n,shape=lam1); zb <- rgamma(n,shape=lam2); 
pm <- exp(0.5+xm*beta+zf)
pf <- exp(0.5+xf*beta+zf)
pf <- pf/(1+pf)
pm <- pm/(1+pm)
pb1 <- exp(0.5+xb1*beta+zf+zb)
pb1 <- pb1/(1+pb1)
pb2 <- exp(0.5+xb2*beta+zf+zb)
pb2 <- pb1/(1+pb2)
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

