plot.cums <-  function (x , pointwise.ci=1, hw.ci=0,
sim.ci=0, robust=0, specific.comps=FALSE,level=0.05, start.time = 0, 
stop.time = 0, add.to.plot=FALSE, mains=TRUE, xlab="Time",
ylab ="Cumulative coefficients",...) 
{
 object<-x; rm(x); 
 B<-object$cum; V<-object$var.cum; p<-dim(B)[[2]]; 
 if (robust>=1) {V<-object$robvar.cum;}

 if (sum(specific.comps)==FALSE) comp<-2:p else comp<-specific.comps+1
 if (stop.time==0) stop.time<-max(B[,1]);

 med<-B[,1]<=stop.time & B[,1]>=start.time
 B<-B[med,]; Bs<-B[1,];  B<-t(t(B)-Bs); B[,1]<-B[,1]+Bs[1];
 V<-V[med,]; Vs<-V[1,]; V<-t( t(V)-Vs); 
 Vrob<-object$robvar.cum; 
 Vrob<-Vrob[med,]; Vrobs<-Vrob[1,]; Vrob<-t( t(Vrob)-Vrobs); 

 c.alpha<- qnorm(1-level/2)
 for (v in comp) { 
 c.alpha<- qnorm(1-level/2)
 est<-B[,v];ul<-B[,v]+c.alpha*V[,v]^.5;nl<-B[,v]-c.alpha*V[,v]^.5;
 if (add.to.plot==FALSE) 
 {
 plot(B[,1],est,ylim=1.05*range(ul,nl),type="s",xlab=xlab,ylab=ylab) 
 if (mains==TRUE) title(main=colnames(B)[v]); }
 else lines(B[,1],est,type="s"); 
 if (pointwise.ci>=1) {
 lines(B[,1],ul,lty=pointwise.ci,type="s");
 lines(B[,1],nl,lty=pointwise.ci,type="s"); }
 if (robust>=1) {
 lines(B[,1],ul,lty=robust,type="s"); 
 lines(B[,1],nl,lty=robust,type="s"); }
 if (hw.ci>=1) {
 if (level!=0.05) cat("Hall-Wellner bands only 95 % \n");
 tau<-length(B[,1])
 nl<-B[,v]-1.13*V[tau,v]^.5*(1+V[,v]/V[tau,v])
 ul<-B[,v]+1.13*V[tau,v]^.5*(1+V[,v]/V[tau,v])
 lines(B[,1],ul,lty=hw.ci,type="s"); 
 lines(B[,1],nl,lty=hw.ci,type="s"); }
 if (sim.ci>=1) {
 if (is.null(object$conf.band)==TRUE) 
 cat("Uniform simulation based bands only computed for n.sim> 0\n")
 if (level!=0.05) c.alpha<-percen(object$sim.testBeq0[,v-1],1-level)
 else c.alpha<-object$conf.band[v-1];
 nl<-B[,v]-c.alpha*Vrob[,v]^.5; ul<-B[,v]+c.alpha*Vrob[,v]^.5;
 lines(B[,1],ul,lty=sim.ci,type="s"); 
 lines(B[,1],nl,lty=sim.ci,type="s"); }
 abline(h=0)
 }
}

plotScore<-function (object,specific.comps=FALSE,mains=TRUE, xlab="Time",ylab ="Cumulative coefficients") 
{
if (ylab=="Cumulative regression function") ylab<-"Test process";

if (class(object)=="cox.aalen") {
            obsProc<-object$test.procProp; 
            simProc<-object$sim.test.procProp; }
else {      obsProc<-object$test.procBeqC; 
            simProc<-object$sim.test.procBeq; }

dim1<-ncol(obsProc)
if (sum(specific.comps)==FALSE) comp<-2:dim1 else comp<-specific.comps+1

for (i in comp) {
ranyl<-range(obsProc[,i]);
for (j in 1:50)
ranyl<-range(c(ranyl,as.matrix(simProc[[j]])[,i-1]));
mr<-max(abs(ranyl)); 

plot(obsProc[,1],obsProc[,i],ylim=c(-mr,mr),lwd=2,xlab=xlab,ylab=ylab,type="s")
if (mains==TRUE) title(main=colnames(obsProc)[i]); 
for (j in 1:50)
lines(obsProc[,1],as.matrix(simProc[[j]])[,i-1],
col="grey",lwd=1,lty=1,type="s")
lines(obsProc[,1],obsProc[,i],lwd=2,type="s")
} 
}
