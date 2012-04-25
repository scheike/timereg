##' @S3method print biprobit
print.biprobit <- function(x,...) {
  printCoefmat(x$coef,...)
  return(invisible(x))
}
