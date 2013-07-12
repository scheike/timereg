##' @S3method plot bptwin
plot.bptwin <- function(x,n=50,rg=range(x$B[,1]),xlab="Time",ylab="Concordance",...) {
  if (x$Blen>0) {
    ##    rg <- range(x$B[,1])
    t <- seq(rg[1],rg[2],length.out=n)
    B0 <- bs(t,degree=x$Blen)
    b0. <- coef(x)[x$midx0]
    b1. <- coef(x)[x$midx1]
    b0 <- trMean(b0.,x$Blen)
    b1 <- trMean(b1.,x$Blen)
    b00 <- tail(b0,x$Blen)
    b11 <- tail(b1,x$Blen)
    pr0 <- sapply(as.numeric(B0%*%b00+b0[1]), function(z)
                  pmvn(upper=rep(z,2),sigma=x$Sigma0))
    pr1 <- sapply(as.numeric(B0%*%b11+b1[1]), function(z)
                  pmvn(upper=rep(z,2),sigma=x$Sigma1))
    plot(pr0~t,type="l", xlab=xlab, ylab=ylab,...)
    lines(pr1~t,type="l",lty=2)
  }
  return(invisible(x))
}
