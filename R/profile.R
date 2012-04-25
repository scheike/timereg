bootstrap <- function (x, ...) UseMethod("bootstrap")
bootstrap.bptwin <- function(x,R=10,silent=FALSE,control=list(),...) {
  mycall <- x$call
  mycall$control <- control
  data <- x$data  
  idx <- seq(nrow(data))
  
  
  bootfun <- function(i,...) {
    if (!silent) cat(".")
    d0 <- data[sample(idx,replace=TRUE),]
    mycall$data=as.name("d0")
    browser()
    e0 <- with(x, eval(mycall))
    coef(e0)[,1]
  }
  if (require(foreach)) {
    res <- foreach (i=seq(R)) %dopar% bootfun(i)    
  }
  else {
    res <- lapply(seq(R),bootfun)
  }
  if (!silent) cat("\n")
  res <- rbind(res,coef(x)[,1])
  return(res)
}


profile.bptwin <- function(fitted,...) {
  mycall <- fitted$call
  mycall$constrain <- c(NA,0.423258,NA)
  mycall$stderr <- FALSE
  mycall$data=as.name("data")
  with(fitted, eval(mycall))
}
