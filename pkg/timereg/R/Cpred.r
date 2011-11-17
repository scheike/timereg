Cpred<-function(cum,xval,Tminus=FALSE,start.val=0,cum.startval=0,order=FALSE) 
{
designX<-as.matrix(cum); 
cumtimes <- designX[,1]
px<-as.integer(dim(designX)[2]);
nx<-as.integer(dim(designX)[1]);
nval<-length(xval); 
pred<-rep(0,nval); 

sout<-.C("Cpred", 
as.double(cumtimes),as.integer(nx),as.integer(px),
as.double(xval),as.integer(nval),as.integer(pred),
as.integer(Tminus),PACKAGE="timereg")

xval.order <- sout[[6]]; 
pred <- predb <-  xval.order
predb <- xval.order
predb[pred==0] <- 1

predcum <- as.matrix(designX[predb,-1])
predcum[pred==0,] <- cum.startval

if (order==FALSE) return(cbind(xval,predcum)) else return(list(xval.order=xval.order,predb=predb))
}
