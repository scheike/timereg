haplomatch.cif<-function(formula,data=sys.parent(),
cause,times,designfuncX,designfuncZ,Nit=50,
clusters=NULL,gamma=0,n.sim=500,weighted=0,model="additive",
causeS=1,cens.code=0,detail=0,interval=0.01,resample.iid=1,
cens.model="KM",time.pow=0,fix.haplofreq=0,haplo.freq=NULL,alpha.iid=NULL,
geno.setup=NULL,fit.haplofreq=NULL,design.test=0,covnamesX=NULL,covnamesZ=NULL)
{
  call <- match.call()
  m <- match.call(expand = FALSE)
out<- haplo.cif( formula = formula,data=data, cause,times,
       designfuncX,designfuncZ, Nit = Nit, clusters=clusters, 
       gamma=gamma, n.sim = n.sim, match=TRUE, 
       weighted=weighted,model=model,causeS=causeS,
       cens.code=cens.code, detail = detail, 
       interval=interval, resample.iid=resample.iid,
       cens.model=cens.model,time.pow=time.pow,
       fix.haplofreq = fix.haplofreq, haplo.freq =haplo.freq,
       alpha.iid=alpha.iid, geno.setup=geno.setup, 
       fit.haplofreq=fit.haplofreq,design.test=design.test,
       covnamesX=covnamesX,covnamesZ=covnamesZ)
  attr(out, "Call") <- sys.call()
  attr(out, "Formula") <- formula
  return(out)
}
