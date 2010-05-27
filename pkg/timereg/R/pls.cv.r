pls.surv.cv<-function(formula=formula(data),data=sys.parent(),
K = 10,start.time=0,max.time=NULL,id=NULL,
pls.dims = seq(from = 0, to = 5, by = 1), 
trace = FALSE, plot.it = FALSE , se = TRUE, silent=1,
    ...) 
{ ## {{{
## {{{ reads design and survival times
    call <- match.call()
    m <- match.call(expand = FALSE)
    m$start.time <- m$max.time <-  m$id <-
    m$K<-m$pls.dims<-m$trace<-m$plot.it<-m$se<- m$silent<- NULL
    special <- c("const")
    Terms <- if (missing(data)) terms(formula, special)
    else terms(formula, special, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    mt <- attr(m, "terms")
    intercept <- attr(mt, "intercept")
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")

   des<-read.design(m,Terms)
   X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px; pz<-des$pz;
   covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ
   clusters <- NULL
   pxz <- px + pz; 
   Yorig<-X[,1]
   X<-X[,-1]

   survs<-read.surv(m,id,npar,clusters,start.time,max.time)
   times<-survs$times;id<-id.call<-survs$id.cal;
   time2<-survs$stop
   status<-survs$status;
## }}}

    antpers = nrow(X); 
    all.folds <- cv.folds(antpers, K)
    residmat <- matrix(0, length(pls.dims), K)

    for (i in seq(K)) {
        omit <- all.folds[[i]]
        Xl<-X[-omit,]; Zl<-Z[-omit,]; 
        datal=data.frame(tt=time2[-omit],status=status[-omit])
        nn<-antpers-length(omit); 
        Xo<-X[omit,]; Zo<-Z[omit,]; 
        datao=data.frame(tt=time2[omit],status=status[omit])
        nno<- length(omit); 
        if (npar==TRUE) 
        outo<-additive.compSs(Surv(tt,status)~const(Xo),data=datao,
            max.time=max.time,start.time=start.time,silent=silent)
        else 
        outo<-additive.compSs(Surv(tt,status)~-1+Xo+const(Zo),data=datao,
            max.time=max.time,start.time=start.time,silent=silent)
        #print(outo$intZHZ); print(outo$intZHdN)
        pfit<-c()
        for (j in pls.dims)
        {
        if (npar==TRUE) 
        out<-additive.pls(Surv(tt,status)~Xl,data=datal,pls.dim=j,
            max.time=max.time,start.time=start.time,silent=silent)
        else out<-additive.pls(Surv(tt,status)~Xl+const(Zl),data=datal,pls.dim=j,
            max.time=max.time,start.time=start.time,silent=silent)
        betas<-matrix(out$tbeta.pls,nrow=1)
        #print(betas)

        fit<- apply(betas,1,function(x)  x %*% outo$intZHZ %*% x) 
        fit<- drop(fit - 2 *betas %*% outo$intZHdN )
        pfit<-c(pfit,fit)
        }
        if (length(omit) == 1) 
            fit <- matrix(pfit, nrow = 1)
        residmat[, i] <- pfit
        if (trace) 
            cat("\n CV Fold", i, "\n\n")
    }
    cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var)/K)
    object <- list(pls.dims = pls.dims, cv = cv, cv.error = cv.error, 
        cv.dim = pls.dims[order(cv)][1])
    if (plot.it) plot(pls.dims,cv,type="l")
    invisible(object)
} ## }}}

cv.folds<- function (n, folds = 10) 
{
    split(sample(1:n), rep(1:folds, length = n))
}

