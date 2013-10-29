##' @export 
sim.bin.fam <- function(n,beta=0.0,theta=1,lam1=1,lam2=1) { ## {{{ 
x1 <- rbinom(n,1,0.5); x2 <- rbinom(n,1,0.5); 
x3 <- rbinom(n,1,0.5); x4 <- rbinom(n,1,0.5); 
###
zf <- rgamma(n,shape=lam1); zb <- rgamma(n,shape=lam2); 
pm <- exp(0.5+x1*beta+zf)
pf <- exp(0.5+x2*beta+zf)
pf <- pf/(1+pf)
pm <- pm/(1+pm)
pb1 <- exp(0.5+x1*beta+zf+zb)
pb1 <- pb1/(1+pb1)
ym <- rbinom(n,1,pm)
yf <- rbinom(n,1,pf)
yb1 <- rbinom(n,1,pb1)
yb2 <- rbinom(n,1,pb1)
#
data.frame(x1=x1,x2=x2,ym=ym,yf=yf,yb1=yb1,yb2=yb2,id=1:n)
} ## }}} 

