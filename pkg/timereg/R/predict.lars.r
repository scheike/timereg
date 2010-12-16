mypredict.lars<-function (object, s, type = c("fit", "coefficients"), 
mode = c("step", "fraction", "norm"), ...) 
{
    newx<-NULL
    mode <- match.arg(mode)
    type <- match.arg(type)
    if (missing(newx) & type == "fit") {
        # warning("Type=fit with no newx argument; type switched to coefficients")
        type <- "coefficients"
    }
    betas <- object$beta
    sbetas <- scale(betas, FALSE, 1/object$normx)
    sbetas<-betas; 
    kp <- dim(betas)
    k <- kp[1]
    p <- kp[2]
    steps <- seq(k)
    if (missing(s)) {
        s <- steps
        mode <- "step"
    }
    sbeta <- switch(mode, step = {
        if (any(s < 0) | any(s > k)) 
            stop("Argument s out of range")
        steps
    }, fraction = {
        if (any(s > 1) | any(s < 0)) 
            stop("Argument s out of range")
        nbeta <- drop(abs(sbetas) %*% rep(1, p))
        nbeta/nbeta[k]
    }, norm = {
        nbeta <- drop(abs(sbetas) %*% rep(1, p))
        if (any(s > nbeta[k]) | any(s < 0)) 
            stop("Argument s out of range")
        nbeta
    })
    sfrac <- (s - sbeta[1])/(sbeta[k] - sbeta[1])
    sbeta <- (sbeta - sbeta[1])/(sbeta[k] - sbeta[1])
    usbeta <- unique(sbeta)
    useq <- match(usbeta, sbeta)
    sbeta <- sbeta[useq]
    betas <- betas[useq, ]
    coord <- approx(sbeta, seq(sbeta), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    newbetas <- ((sbeta[right] - sfrac) * betas[left, , drop = FALSE] + 
        (sfrac - sbeta[left]) * betas[right, , drop = FALSE])/(sbeta[right] - 
        sbeta[left])
    newbetas[left == right, ] <- betas[left[left == right], ]
    fitf<- apply(newbetas,1,function(x)  x %*% object$Gram %*% x) 
    fitf<- fitf - 2 *newbetas %*% object$intZHdN 
    robject <- switch(type, coefficients = list(s = s, fraction = sfrac, 
                mode = mode, coefficients = drop(newbetas)), 
                fit = list(s = s, fraction = sfrac, mode = mode, fit = fitf)) 
robject
}

mycoef.lars<- function (object, ...)
{
    mypredict.lars(object, type = "coefficient", ...)$coef
}

myplot.lars<-
function (x, xvar = c("norm", "df", "arc.length"), breaks = TRUE, 
    plottype = c("coefficients", "Cp"), omit.zeros = TRUE, eps = 1e-10, 
    ...) 
{
    object <- x
    plottype <- match.arg(plottype)
    xvar <- match.arg(xvar)
    coef1 <- object$beta
    # coef1 <- scale(coef1, FALSE, 1/object$normx)
    if (omit.zeros) {
        c1 <- drop(rep(1, nrow(coef1)) %*% abs(coef1))
        nonzeros <- c1 > eps
        cnums <- seq(nonzeros)[nonzeros]
        coef1 <- coef1[, nonzeros]
    }
    else cnums <- seq(ncol(coef1))
    s1 <- switch(xvar, norm = {
        s1 <- apply(abs(coef1), 1, sum)
        s1/max(s1)
    }, df = seq(length(object$arc.length) + 1), arc.length = cumsum(c(0, 
        object$arc.length)))
    xname <- switch(xvar, norm = "|beta|/max|beta|", df = "Df", 
        arc.length = "Arc Length")
    if (plottype == "Cp") {
        Cp <- object$Cp
        plot(s1, Cp, type = "b", xlab = xname, main = object$type, 
            ...)
    }
    else {
        matplot(s1, coef1, xlab = xname, ..., type = "b", pch = "*", 
            ylab = "Coefficients")
        title(object$type, line = 2.5)
        abline(h = 0, lty = 3)
        axis(4, at = coef1[nrow(coef1), ], labels = paste(cnums), 
            cex = 0.8, adj = 0)
        if (breaks) {
            axis(3, at = s1, labels = paste(seq(s1) - 1), cex = 0.8)
            abline(v = s1)
        }
    }
    invisible()
}
