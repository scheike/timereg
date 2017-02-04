##' Set global options for \code{mets}
##'
##' Extract and set global parameters of \code{mets}. 
##'
##' \itemize{
##'   \item \code{regex}: If TRUE character vectors will be interpreted as regular expressions (\code{dby}, \code{dcut}, ...)
##'   \item \code{silent}: Set to \code{FALSE} to disable various output messages
##' }
##' @param ... Arguments
##' @return \code{list} of parameters
##' @keywords models
##' @examples
##' \dontrun{
##' mets.options(regex=TRUE)
##' }
##' @export
mets.options <- function(...) {
    dots <- list(...)
    newopt <- curopt <- get("options",envir=mets.env)
    if (length(dots)==0)
        return(curopt)
    if (length(dots)==1 && is.list(dots[[1]]) && is.null(names(dots))) {
        dots <- dots[[1]]
    }
    idx <- which(names(dots)!="")
    newopt[names(dots)[idx]] <- dots[idx]
    assign("options",newopt,envir=mets.env)
    invisible(curopt)
}

mets.env <- new.env()
assign("options",
       list(debug=FALSE,
            regex=FALSE
            ), envir=mets.env)
