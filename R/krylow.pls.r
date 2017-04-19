###krylow.pls<-function(D,d,dim)
###{
###R=d;  Sxxsxy=R;
###if (dim>=2)
###for (i in 2:dim)
###{
###Sxxsxy=D %*% Sxxsxy  ;
###R=cbind(R,Sxxsxy);
###}
###beta= R %*% solve(t(R) %*% D %*% R) %*% t(R) %*% d
###beta=  beta;
###return(list(beta=beta))
###}

### Anders Gorst Rasmussen's code


#' Fits Krylow based PLS for additive hazards model
#' 
#' Fits the PLS estimator for the additive risk model based on the least
#' squares fitting criterion
#' 
#' \deqn{ L(\beta,D,d) = \beta^T D \beta - 2 \beta^T d } where \eqn{D=\int Z H
#' Z dt} and \eqn{d=\int Z H dN}.
#' 
#' 
#' @param D defined above
#' @param d defined above
#' @param dim number of pls dimensions
#' @return returns a list with the following arguments: \item{beta}{PLS
#' regression coefficients}
#' @author Thomas Scheike
#' @references Martinussen and Scheike, The Aalen additive hazards model with
#' high-dimensional regressors, submitted.
#' 
#' Martinussen and Scheike, Dynamic Regression Models for Survival Data,
#' Springer (2006).
#' @keywords survival
#' @examples
#' 
#' ## makes data for pbc complete case
#' data(mypbc)
#' pbc<-mypbc
#' pbc$time<-pbc$time+runif(418)*0.1; pbc$time<-pbc$time/365
#' pbc<-subset(pbc,complete.cases(pbc));
#' covs<-as.matrix(pbc[,-c(1:3,6)])
#' covs<-cbind(covs[,c(1:6,16)],log(covs[,7:15]))
#' 
#' ## computes the matrices needed for the least squares 
#' ## criterion 
#' out<-aalen(Surv(time,status>=1)~const(covs),pbc,robust=0,n.sim=0)
#' S=out$intZHZ; s=out$intZHdN;
#' 
#' out<-krylow.pls(S,s,dim=2)
#' 
#' @export
krylow.pls <- function(D,d,dim = 1) {
	r <- d
	p <- r
	rsold <- drop(t(r) %*% r) 
	x <- rep(0,length(d)) 
	for(k in 1:dim){
		Ap <- D %*% p;
		alpha <-drop(rsold / (t(p) %*% Ap)); x <- x + alpha * p;
		r <- r - alpha * Ap;
		rsnew <- drop(t(r) %*% r);
		p <- r + rsnew / rsold * p;
		rsold <-rsnew;
	}
	return(x) 
}
