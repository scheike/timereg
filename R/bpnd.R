bpnd <- function(u1,u2,r) {
  arglist <- list(name="cdf",
                  a=u1,
                  b=u2,
                  r=r,
                  DUP=FALSE,PACKAGE="bptwin")
  res <- do.call(".Call",arglist)
  return(res)
}
