#' Cut a variable
#' 
#' Calls the cut function to cut variables on data frame.
#' 
#' 
#' @param x variable to cut
#' @param cuts number of groups, 4 gives quartiles
#' @param breaks can also give breaks
#' @param ...  other argument for cut function of R
#' @author Thomas Scheike
#' @keywords survival
#' @examples
#' 
#' data(sTRACE)
#' gx <- qcut(sTRACE$age)
#' table(gx)
#' 
#' @export
qcut <- function(x,cuts=4,breaks=NULL,...)
{# {{{
	if (is.null(breaks)) {
		   probs <- seq(0,1,length.out=cuts+1)
	   bb <- quantile(x,probs)
	} else bb <- breaks
	gx<- cut(x,breaks=bb,include.lowest=TRUE,...)
	return(gx)
}# }}}

