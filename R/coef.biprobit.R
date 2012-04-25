##' @S3method coef biprobit
coef.biprobit <- function(object,matrix=FALSE,...) {
  if (matrix) return(object$coef)
  return(object$coef[,1])
}
