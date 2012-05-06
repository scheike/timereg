##' Analysis of Multivariate Events 
##' 
##' Implementation of various statistical models for multivariate
##' event history data. Including multivariate cumulative incidence models,
##' and bivariate random effects probit models (Liability models) 
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

##' np data set
##'
##' @name np
##' @docType data
##' @references \url{data_blah.com}
##' @keywords data
NULL

##' multcif data set
##'
##' @name multcif
##' @docType data
##' @references \url{data_blah.com}
##' @keywords data
NULL

##' Method for extracting score of model object
##'
##' @title Method for extracting score of model object
##' @param x Model object
##' @param ... Additional arguments
##' @author Klaus K. Holst
##' @export
score <- function(x,...) UseMethod("score")
