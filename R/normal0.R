
##' @export
loglikMVN <- function(yl,yu,status,mu,S,thres) {
    .Call("loglikMVN",
          yl=as.matrix(yl),
          yu=as.matrix(yu),
          status=as.integer(status),
          mu=as.matrix(mu),dmu=NULL,s=as.matrix(S),ds=NULL,
          z=NULL,su=NULL,dsu=NULL,
          threshold=as.matrix(thres),
          dthreshold=NULL, package="mets")        
}
