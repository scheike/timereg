
#' @export
residualsTimereg <- function(object,data=data)
{ ## {{{
### computes residuals for data based on model given in object 

if (!inherits(object,c("cox.aalen","aalen"))) stop("Computes residuals for Aalen or Cox.aalen object") 
else {
 formula<-attr(object,"Formula");
 beta.fixed <- attr(object,"beta.fixed")
 if (is.null(beta.fixed)) beta.fixed <- 1; 
 model <- class(object); 
 ldata<-aalen.des(formula,data=data,model=model);
 id <- attr(object,"id"); 
 mclusters <- attr(object,"cluster")
 X<-ldata$X; 
 time2<-ldata$time2; 
 start<-ldata$time; 
 Z<-ldata$Z;  
 status<-ldata$status;
 otime2 <- attr(object,"stop"); 
 ostart <- attr(object,"start");
 ostatus <- attr(object,"status");

 if (!is.null(attr(object,"max.time"))) status <- status*(time2< attr(object,"max.time")); 
 antpers<-nrow(X);
 if (is.null(Z)==TRUE) {npar<-TRUE; semi<-0;}  else { Z<-as.matrix(Z); npar<-FALSE; semi<-1;}
 if (npar==TRUE) {Z<-matrix(0,antpers,1); pz<-1; fixed<-0;} else {fixed<-1;pz<-ncol(Z);}
 px<-ncol(X);

 if (sum(abs(start))>0) lefttrunk <- 1  else lefttrunk <- 0;  
 cumhazleft <- 0; 
 nn <- nrow(object$cum) 

 cum <- Cpred(object$cum,time2)[,-1]
 cumhaz0 <- apply(cum*X,1,sum)
 cumhazleft <- rep(0,antpers)
 RR <- rep(1,antpers); 

if (inherits(object,"cox.aalen"))
{ ## {{{
  RR <- exp(Z %*% object$gamma); 
  cumhaz <- cumhaz0 * RR;
  if (lefttrunk==1) {
      cum <- Cpred(object$cum,start)[,-1]
      cumhazleft <- apply(cum*X,1,sum)
      cumhazleft <- cumhazleft * RR;
  }
} ## }}}

if (inherits(object,"aalen"))
{#{{{
  if (npar==FALSE) { ## semi-parametric risk model
      ex.haz <- (Z %*% object$gamma) ; 
      cumhaz <- cumhaz0+ex.haz*time2
     if (lefttrunk==1) {
	 cum <- Cpred(object$cum,start)[,-1]
	 cumhazleft <- apply(cum*X,1,sum)
	 cumhazleft  <-  cumhazleft+ex.haz*start
     }
  } else {  ## Aalen model
	  cumhaz <- cumhaz0
          if (lefttrunk==1) {
	     cum <- Cpred(object$cum,start)[,-1]
	     cumhazleft <- apply(cum*X,1,sum)
	     if (npar==TRUE) cumhazleft <-  cumhazleft
          }
  }
} #}}}

} 

residuals <- status- cumhaz
out <- list(residuals=c(residuals),status=c(status),cumhaz=c(cumhaz),cumhazleft=c(cumhazleft),RR=RR)
} ## }}}

