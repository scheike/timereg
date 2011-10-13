Cpred <- function(cum,xval){
  designX <- as.matrix(cum)
  px <- as.integer(dim(designX)[2])
  nx <- as.integer(dim(designX)[1])
  nval <- length(xval)
  pred <- matrix(0,nval,px)
  sout <- .C("Cpred",
             as.double(cum),
             as.integer(nx),
             as.integer(px),
             as.double(xval),
             as.integer(nval),
             pred=as.double(pred),
             PACKAGE="compRisk")
  pred <- matrix(sout$pred,nval,px)
  return(pred)
}
