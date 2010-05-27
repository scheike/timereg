my.lars<-function (Gram, xy, n, y, type = c("lasso", "lar", "forward.stagewise", 
    "stepwise"), trace = FALSE, normalize = FALSE, intercept = FALSE , 
     eps = .Machine$double.eps, max.steps, use.Gram = TRUE) 
{
    if (!require("lars")) stop("package Lars not available")
    y<-rep(1,n); x<-rep(1,n)
    call <- match.call()
    type <- match.arg(type)
    TYPE <- switch(type, lasso = "LASSO", lar = "LAR", forward.stagewise = "Forward Stagewise", 
        stepwise = "Forward Stepwise")
    if (trace) 
        cat(paste(TYPE, "sequence\n"))
    #nm <- dim(x) #n <- nm[1] #m <- nm[2]
    m <- ncol(Gram)
    im <- inactive <- seq(m)
    one <- rep(1, n)
    vn <- dimnames(x)[[2]]
    if (intercept) {
        meanx <- drop(one %*% x)/n
        x <- scale(x, meanx, FALSE)
        mu <- mean(y)
        y <- drop(y - mu)
    }
    else {
        meanx <- rep(0, m)
        mu <- 0
        y <- drop(y)
    }
    if (normalize) {
        normx <- sqrt(drop(one %*% (x^2)))
        nosignal <- normx/sqrt(n) < eps
        if (any(nosignal)) {
            ignores <- im[nosignal]
            inactive <- im[-ignores]
            normx[nosignal] <- eps * sqrt(n)
            if (trace) 
                cat("LARS Step 0 :\t", sum(nosignal), "Variables with Variance < eps; dropped for good\n")
        }
        else ignores <- NULL
        names(normx) <- NULL
        x <- scale(x, FALSE, normx)
    }
    else {
        normx <- rep(1, m)
        ignores <- NULL
    }
    if (use.Gram & missing(Gram)) {
        if (m > 500 && n < m) 
            cat("There are more than 500 variables and n<m;\nYou may wish to restart and set use.Gram=FALSE\n")
        if (trace) 
            cat("Computing X'X .....\n")
        Gram <- t(x) %*% x
    }
    Cvec<-xy;
    ssy <- sum(y^2)
    residuals <- y
    if (missing(max.steps)) 
        max.steps <- 8 * min(m, n - intercept)
    beta <- matrix(0, max.steps + 1, m)
    lambda = double(max.steps)
    Gamrat <- NULL
    arc.length <- NULL
    R2 <- 1
    RSS <- ssy
    first.in <- integer(m)
    active <- NULL
    actions <- as.list(seq(max.steps))
    drops <- FALSE
    Sign <- NULL
    R <- NULL
    k <- 0
    while ((k < max.steps) & (length(active) < min(m - length(ignores), 
        n - intercept))) {
        action <- NULL
        C <- Cvec[inactive]
        Cmax <- max(abs(C))
        if (Cmax < eps * 100) {
            if (trace) 
                cat("Max |corr| = 0; exiting...\n")
            break
        }
        k <- k + 1
        lambda[k] = Cmax
        if (!any(drops)) {
            new <- abs(C) >= Cmax - eps
            C <- C[!new]
            new <- inactive[new]
            for (inew in new) {
                if (use.Gram) {
                  R <- lars::updateR(Gram[inew, inew], R, drop(Gram[inew, 
                    active]), Gram = TRUE, eps = eps)
                }
                else {
                  R <- lars::updateR(x[, inew], R, x[, active], Gram = FALSE, 
                    eps = eps)
                }
                if (attr(R, "rank") == length(active)) {
                  nR <- seq(length(active))
                  R <- R[nR, nR, drop = FALSE]
                  attr(R, "rank") <- length(active)
                  ignores <- c(ignores, inew)
                  action <- c(action, -inew)
                  if (trace) 
                    cat("LARS Step", k, ":\t Variable", inew, 
                      "\tcollinear; dropped for good\n")
                }
                else {
                  if (first.in[inew] == 0) 
                    first.in[inew] <- k
                  active <- c(active, inew)
                  Sign <- c(Sign, sign(Cvec[inew]))
                  action <- c(action, inew)
                  if (trace) 
                    cat("LARS Step", k, ":\t Variable", inew, 
                      "\tadded\n")
                }
            }
        }
        else action <- -dropid
        Gi1 <- backsolve(R, lars::backsolvet(R, Sign))
        dropouts <- NULL
        if (type == "forward.stagewise") {
            directions <- Gi1 * Sign
            if (!all(directions > 0)) {
                if (use.Gram) {
                  nnls.object <- lars::nnls.lars(active, Sign, R, directions, 
                    Gram[active, active], trace = trace, use.Gram = TRUE, 
                    eps = eps)
                }
                else {
                  nnls.object <- lars::nnls.lars(active, Sign, R, directions, 
                    x[, active], trace = trace, use.Gram = FALSE, 
                    eps = eps)
                }
                positive <- nnls.object$positive
                dropouts <- active[-positive]
                action <- c(action, -dropouts)
                active <- nnls.object$active
                Sign <- Sign[positive]
                Gi1 <- nnls.object$beta[positive] * Sign
                R <- nnls.object$R
                C <- Cvec[-c(active, ignores)]
            }
        }
        A <- 1/sqrt(sum(Gi1 * Sign))
        w <- A * Gi1
        if (!use.Gram) 
            u <- drop(x[, active, drop = FALSE] %*% w)
        if ((length(active) >= min(n - intercept, m - length(ignores))) | 
            type == "stepwise") {
            gamhat <- Cmax/A
        }
        else {
            if (use.Gram) {
                a <- drop(w %*% Gram[active, -c(active, ignores), 
                  drop = FALSE])
            }
            else {
                a <- drop(u %*% x[, -c(active, ignores), drop = FALSE])
            }
            gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))
            gamhat <- min(gam[gam > eps], Cmax/A)
        }
        if (type == "lasso") {
            dropid <- NULL
            b1 <- beta[k, active]
            z1 <- -b1/w
            zmin <- min(z1[z1 > eps], gamhat)
            if (zmin < gamhat) {
                gamhat <- zmin
                drops <- z1 == zmin
            }
            else drops <- FALSE
        }
        beta[k + 1, ] <- beta[k, ]
        beta[k + 1, active] <- beta[k + 1, active] + gamhat * 
            w
        if (use.Gram) {
            Cvec <- Cvec - gamhat * Gram[, active, drop = FALSE] %*% 
                w
        }
        else {
            residuals <- residuals - gamhat * u
            Cvec <- drop(t(residuals) %*% x)
        }
        Gamrat <- c(Gamrat, gamhat/(Cmax/A))
        arc.length <- c(arc.length, gamhat)
        if (type == "lasso" && any(drops)) {
            dropid <- seq(drops)[drops]
            for (id in rev(dropid)) {
                if (trace) 
                  cat("Lasso Step", k + 1, ":\t Variable", active[id], 
                    "\tdropped\n")
                R <- lars::downdateR(R, id)
            }
            dropid <- active[drops]
            beta[k + 1, dropid] <- 0
            active <- active[!drops]
            Sign <- Sign[!drops]
        }
        if (!is.null(vn)) 
            names(action) <- vn[abs(action)]
        actions[[k]] <- action
        inactive <- im[-c(active, ignores)]
        if (type == "stepwise") 
            Sign = Sign * 0
    }
    beta <- beta[seq(k + 1), , drop = FALSE]
    lambda = lambda[seq(k)]
    dimnames(beta) <- list(paste(0:k), vn)
    if (trace) 
        cat("Computing residuals, RSS etc .....\n")
    #residuals <- y - x %*% t(beta)
    residuals <- 0 
    beta <- scale(beta, FALSE, normx)
    RSS <- 0 # apply(residuals^2, 2, sum)
    R2 <- 0 #1 - RSS/RSS[1]
    actions = actions[seq(k)]
    netdf = sapply(actions, function(x) sum(sign(x)))
    df = cumsum(netdf)
    if (intercept) 
        df = c(Intercept = 1, df + 1)
    else df = c(Null = 0, df)
    rss.big = rev(RSS)[1]
    df.big = n - rev(df)[1]
    if (rss.big < eps | df.big < eps) 
        sigma2 = NaN
    else sigma2 = rss.big/df.big
    Cp <- RSS/sigma2 - n + 2 * df
    attr(Cp, "sigma2") = sigma2
    attr(Cp, "n") = n
    object <- list(call = call, type = TYPE, df = df, lambda = lambda, 
        R2 = R2, RSS = RSS, Cp = Cp, actions = actions[seq(k)], 
        entry = first.in, Gamrat = Gamrat, arc.length = arc.length, 
        Gram = if (use.Gram) Gram else NULL, beta = beta, mu = mu, 
        intZHdN=xy, normx = normx, meanx = meanx)
    class(object) <- "lars"
    object
}
