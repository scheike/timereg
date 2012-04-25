##' Analysis of Multivariate Events 
##' 
##' Implementation of various statistical models for multivariate
##' event history data. Including multivariate cumulative incidence models,
##' and bivariate random effects probit models (Liability models) 
##'
##' 
##' @name mets-package
##' @docType package
##' @author Klaus K. Holst and Thomas Scheike
##' @useDynLib mets
##' @import stats
##' @keywords package
##' @examples
##' 
##' ## To appear
##' 
NULL

##' @export
`score` <- function(x,...) UseMethod("score")
