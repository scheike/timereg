Cpred<-function(cum,xval) 
{
designX<-as.matrix(cum); 
px<-as.integer(dim(designX)[2]);
nx<-as.integer(dim(designX)[1]);
nval<-length(xval); 
pred<-matrix(0,nval,px); 

sout<-.C("Cpred", 
as.double(cum),as.integer(nx),as.integer(px),
as.double(xval),as.integer(nval),as.double(pred),PACKAGE="timereg")

pred<-matrix(sout[[6]],nval,px);
return(pred)
}
