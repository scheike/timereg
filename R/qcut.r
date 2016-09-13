
qcut <- function(x,cuts=4,...)
{
probs <- seq(0,1,length.out=cuts+1)
bb <- quantile(x,probs,...)
gx<- cut(x,breaks=bb,include.lowest=TRUE)
return(gx)
}


