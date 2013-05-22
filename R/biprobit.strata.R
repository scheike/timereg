do.twinlm.strata <- function(x,fun,...) {
  res <- lapply(x$model,function(m) do.call(fun,c(list(m),list(...))))
  names(res) <- names(x$model)
  class(res) <- "do.twinlm.strata"
  return(res)
}

##' @S3method print do.twinlm.strata
print.do.twinlm.strata <- function(x,...) {
  for (i in seq_len(length(x))) {    
    message(rep("-",60),sep="")
    message("Strata '",names(x)[i],"'",sep="")
    print(x[[i]])
  }
  return(invisible(x))
}


##' @S3method plot twinlm.strata
plot.twinlm.strata <- function(x,...)
  suppressMessages(do.twinlm.strata(x,"plot",...))

##' @S3method print twinlm.strata
print.twinlm.strata <- function(x,...)
  print.do.twinlm.strata(x$model,...)

##' @S3method summary twinlm.strata
summary.twinlm.strata <- function(object,...)
  do.twinlm.strata(object,"summary",...)

##' @S3method coef twinlm.strata
coef.twinlm.strata <- function(object,...) object$coef

##' @S3method logLik twinlm.strata
logLik.twinlm.strata <- function(object,indiv=FALSE,list=FALSE,...) {
  ll <- lapply(object$model,function(x) logLik(x,indiv=indiv,...))
  if (list) return(ll)
  if (!indiv) {
    res <- structure(sum(unlist(ll)),df=0,nall=0)
    for (i in seq(length(ll))) {
      attributes(res)$nall <- attributes(res)$nall+attributes(ll[[i]])$nall
      attributes(res)$df <- attributes(res)$df+attributes(ll[[i]])$df
    }
    ##  attributes(res)$nobs <- attributes(res)$nall-attributes(res)$df
    attributes(res)$nobs <- attributes(res)$nall
    class(res) <- "logLik"
    return(res)
  }
  return(unlist(ll))
}

##' @S3method score twinlm.strata
score.twinlm.strata <- function(x,...) {
  ss <- lapply(x$model,function(m) score(m,indiv=FALSE,...))
  return(unlist(ss))
}
