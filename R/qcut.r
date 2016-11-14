
qcut <- function(x,cuts=4,breaks=NULL,...)
{# {{{
	if (is.null(breaks)) {
		   probs <- seq(0,1,length.out=cuts+1)
	   bb <- quantile(x,probs,...)
	} else bb <- breaks
	gx<- cut(x,breaks=bb,include.lowest=TRUE)
	return(gx)
}# }}}

