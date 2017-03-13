haplo.match<- function ( formula = formula(data),data=sys.parent(),
    designfuncX,designfuncZ, beta = 0,  
    Nit = 10, detail = 0, start.time = 0, max.time = NULL, id = NULL, n.sim = 500, 
    geno.type = NULL, geno.setup=NULL, haplo.freq = NULL,
    fix.beta = 0, fix.haplofreq = 0, two.stage = 0, weighted.test = 0,
    step = 1, lev.marq=1, min.lev.marq=0,
    haplo.design=NULL,haplo.baseline=NULL,alpha=NULL,resample.iid=1,
    covnamesX=NULL,covnamesZ=NULL)
{
  call <- match.call()
  m <- match.call(expand = FALSE)
  out<-haplo.surv( formula = formula,data=data,
    designfuncX,designfuncZ, beta = beta,  match=TRUE, 
    Nit = Nit, detail = detail, start.time = start.time, 
    max.time = max.time, id =id, n.sim = n.sim, 
    geno.type = geno.type, geno.setup=geno.setup, haplo.freq =haplo.freq,
    fix.beta = fix.beta, fix.haplofreq = fix.haplofreq, 
    two.stage = two.stage, weighted.test = weighted.test,
    step = step, lev.marq=lev.marq, min.lev.marq=min.lev.marq,
    haplo.design=haplo.design,haplo.baseline=haplo.baseline,
    alpha=alpha,resample.iid=resample.iid,
    covnamesX=covnamesX,covnamesZ=covnamesZ)
  attr(out, "Call") <- sys.call()
  attr(out, "Formula") <- formula
  return(out)

}

