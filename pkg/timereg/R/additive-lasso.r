surv.lars<-function(S,s,n, l1.weights=NULL, ...) 
{
    if (!require("lars")) stop("package Lars not available")
    if (is.null(l1.weights)==FALSE) adap<-1 else adap<-0; 
    adap<-0; 
    if (adap==1) S<-t(t(S)/l1.weights)
    fit <- my.lars(S, s ,n , ...)
    if (adap==1) fit$coef<-t(t(fit$coef)/l1.weights)
return(fit)
}

surv.lars.cv<-function(formula=formula(data),data=sys.parent(),
l1.weights = NULL, K = 10,start.time=0,max.time=NULL,id=NULL,
fraction = seq(from = 0, to = 1, length = 100), 
trace = FALSE, plot.it = FALSE , se = TRUE, silent=1, ...) 
{
   if (!require("lars")) stop("package Lars not available")
## {{{ reads design and survival times
    call <- match.call()
    m <- match.call(expand = FALSE)
    m$start.time <- m$max.time <-  m$id <-m$silent <- 
    m$covs<-m$l1.weight<-m$K<-m$fraction<-m$trace<-m$plot.it<-m$se<- NULL
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

   survs<-read.surv(m,id,npar,clusters,start.time,max.time)
   times<-survs$times;id<-id.call<-survs$id.cal;
   time2<-survs$stop
   status<-survs$status;
## }}}

    antpers = nrow(Z); 
    all.folds <- cv.folds(antpers, K)
    residmat <- matrix(0, length(fraction), K)
    if (is.null(l1.weights) == FALSE) 
        adap <- 1
    else adap <- 0
    for (i in seq(K)) {
        omit <- all.folds[[i]]
        Xl<-X[-omit,]; Zl<-Z[-omit,]; 
        datal=data.frame(tt=time2[-omit],status=status[-omit])
        nn<-antpers-length(omit); 
        out<-additive.compSs(Surv(tt,status)~-1+Xl+const(Zl),data=datal,
            max.time=max.time,start.time=start.time,silent=silent)
        ###S=out$intZHZ; s=out$intZHdN
        fitno<-my.lars(out$intZHZ,out$intZHdN,nn)
        betas<-mypredict.lars(fitno, fraction, type="coefficients",
               mode = "fraction")$coef

        Xo<-X[omit,]; Zo<-Z[omit,]; 
        datao=data.frame(tt=time2[omit],status=status[omit])
        nno<- length(omit); 
        outo<-additive.compSs(Surv(tt,status)~-1+Xo+const(Zo),data=datao,
            max.time=max.time,start.time=start.time,silent=silent)

        fit<- apply(betas,1,function(x)  x %*% outo$intZHZ %*% x) 
        fit<- drop(fit - 2 *betas %*% outo$intZHdN )

        if (length(omit) == 1) 
            fit <- matrix(fit, nrow = 1)
        residmat[, i] <- fit
        if (trace) 
            cat("\n CV Fold", i, "\n\n")
    }
    cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var)/K)
    object <- list(fraction = fraction, cv = cv, cv.error = cv.error, 
        cv.frac = fraction[order(cv)][1])
    if (plot.it) 
        lars::plotCVLars(object, se = se)
    invisible(object)
}
