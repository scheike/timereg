##' @title Fast reshape
##' Simple reshape/tranpose of data
##' @param data data.frame or matrix
##' @param id id-variable. If omitted then reshape Wide->Long. 
##' @param varying Vector of prefix-names of the time varying variables. Optional for Long->Wide reshaping.
##' @param num Optional number/time variable
##' @param sep String seperating prefix-name with number/time
##' @param keep Vector of column names to keep
##' @param idname Name of id-variable (Wide->Long)
##' @param numname Name of number-variable (Wide->Long)
##' @param factors.keep If false all factors are converted to integers
##' @param idcombine If TRUE and \code{id} is vector of several variables, the unique id is combined from all the variables.
##' Otherwise the first variable is only used as identifier.
##' @param labelnum If TRUE varying variables in wide format (going from long->wide) are labeled 1,2,3,... otherwise use 'num' variable. In long-format (going from wide->long) varying variables matching 'varying' prefix are only selected if their postfix is a number.
##' @param ... Optional additional arguments
##' @author Thomas Scheike, Klaus K. Holst
##' @export
##' @examples
##' library(lava)
##' m <- lvm(c(y1,y2,y3,y4)~x)
##' d <- sim(m,5)
##' fast.reshape(fast.reshape(d,var="y"),id="id")
##' 
##' ##### From wide-format
##' (dd <- fast.reshape(d,var="y"))
##' ## Same with explicit setting new id and number variable/column names
##' ## and seperator "" (default) and dropping x
##' fast.reshape(d,var="y",idname="a",timevar="b",sep="",keep=c())
##' ## Same with 'reshape' list-syntax
##' fast.reshape(d,var=list(c("y1","y2","y3","y4")))
##' 
##' ##### From long-format
##' fast.reshape(dd,id="id")
##' ## Restrict set up within-cluster varying variables
##' fast.reshape(dd,id="id",var="y")
##' fast.reshape(dd,id="id",var="y",keep="x",sep=".")
##' 
##' #####
##' x <- data.frame(id=c(5,5,6,6,7),y=1:5,x=1:5,tv=c(1,2,2,1,2))
##' (xw <- fast.reshape(x,id="id"))
##' (xl <- fast.reshape(xw,varying=c("y","x"),idname="id2",keep=c()))
##' (xl <- fast.reshape(xw,varying=c("y","x","tv")))
##' (xw2 <- fast.reshape(xl,id="id",num="num"))
##' fast.reshape(xw2,varying=c("y","x"),idname="id")
##' 
##' ### more generally:
##' ### varying=list(c("ym","yf","yb1","yb2"), c("zm","zf","zb1","zb2"))
##' ### varying=list(c("ym","yf","yb1","yb2")))
##' 
##' ##### Family cluster example
##' d <- mets:::simBinFam(3)
##' d
##' dd <- fast.reshape(d,var="y")
##' dd
##' dd2 <- fast.reshape(d,varying=list(c("ym","yf","yb1","yb2")))
##' dd2
##' ##'
##' ##'
##' 
##' d <- sim(lvm(~y1+y2+ya),10)
##' (dd <- fast.reshape(d,varying="y"))
##' fast.reshape(d,varying="y",labelnum=TRUE)
##' fast.reshape(dd,id="id",num="num")
##' fast.reshape(dd,id="id",num="num",labelnum=TRUE)
##' fast.reshape(d,vary=c(a="y"),labelnum=TRUE) ## New column name
##' 
##' ##### Prostate cancer example
##' data(prt)
##' head(prtw <- fast.reshape(prt,"id",var="cancer"))
##' ftable(cancer1~cancer2,data=prtw)
##' rm(prtw)
fast.reshape <- function(data,varying,id,num,sep="",keep,
                         idname="id",numname="num",factors.keep=TRUE,
                         idcombine=TRUE,labelnum=FALSE,...) {
    if (!is.data.frame(data) & is.list(data)) {
        data <- as.data.frame(data)
    } else {
        if (NCOL(data)==1) data <- cbind(data)
    }
    if (missing(id)) {
        ## reshape from wide to long format. 
        nn <- colnames(data)
        nsep <- nchar(sep)
        if (missing(varying)) stop("Prefix of time-varying variables needed")
        varsubst <- substitute(varying)
        if (as.character(varsubst)[1]=="-") {
            notvarying <- varsubst[[-1]]
            vars <- setdiff(nn,eval(notvarying))
            browser()
        }
        ## nl <- as.list(seq_along(data)); names(nl) <- nn
        ## varying <- eval(substitute(varying),nl,parent.frame())
        vnames <- NULL
        ncvar <- sapply(varying,nchar)
        newlist <- c()
        numlev <- TRUE
        all_levels <- c()
        thelevels <- c()
        if (!is.list(varying)) {
            for (i in seq_len(length(varying))) {
                ii <- which(varying[i]==substr(nn,1,ncvar[i]))                
                thelevel <- substring(nn[ii],ncvar[i]+1+nsep)
                all_levels <- union(all_levels,thelevel)
                if (labelnum) {
                    ii0 <- suppressWarnings(which(!is.na(as.numeric(thelevel))))
                    ii <- ii[ii0]
                    thelevel <- thelevel[ii0]
                }
                thelevels <- c(thelevels,list(thelevel))
                suppressWarnings(tt <- as.numeric(thelevel))
                newlist <- c(newlist,list(nn[ii[order(tt)]]))
            }
            len <- unlist(lapply(newlist,length))
            for (i in seq_len(length(varying))) {
                if (len[i]<length(all_levels)) {
                    pp <- setdiff(all_levels,thelevels[[i]])
                    vv <- paste(varying[i],pp,sep=sep)
                    data[,vv] <- NA
                    newlist[[i]] <- paste(varying[i],all_levels,sep=sep)
                }                    
            }
            thelevels <- all_levels
            if (any(is.na(as.numeric(all_levels)))) {
                numlev <- FALSE
            }            
            vnames0 <- names(varying) 
            vnames <- varying
            if (!is.null(vnames0)) {
                vidx <- which(vnames0!="")
                vnames[vidx] <- vnames0[vidx]
            }
            varying <- newlist
        } else {
            thelevels <- seq(length(varying[[1]]))
        }
        is_df <- is.data.frame(data)
        oldreshape <- FALSE
        if (is_df) {
            ##      D0 <- droplevels(data)[1,,drop=FALSE]
            D0 <- data[1,,drop=FALSE]
            classes <- unlist(lapply(D0,class))
            dim <- unlist(lapply(D0,NCOL))
            if (any(dim>1) || !all(classes%in%c("numeric","logical","factor","integer","matrix"))) { ## e.g. Surv columns 
                oldreshape <- TRUE
            } else {
                data <- data.matrix(data)
            }
        }
        if (is.null(vnames)) {
            vnames <- unlist(lapply(varying,function(x) x[1]))
            if (!is.null(names(vnames))) vnames <- names(vnames)
        }

        if (oldreshape) return(reshape(as.data.frame(data),varying=varying,direction="long",v.names=vnames,...)) ### Fall-back to stats::reshape

        fixed <- setdiff(nn,unlist(c(varying,numname)))
        if (!missing(keep)) fixed <- intersect(fixed,c(keep,idname,numname))
        nfixed <- length(fixed)
        nvarying <- length(varying)
        nclusts <- unlist(lapply(varying,length))
        ##        nclust <- length(varying[[1]])
        nclust <- max(nclusts)
        if (any(nclusts!=nclust)) stop("Different length of varying vectors!")
        data <- data[,c(fixed,unlist(varying)),drop=FALSE]
        
        long <- as.data.frame(.Call("FastLong",
                                    idata=data,
                                    inclust=as.integer(nclust),
                                    as.integer(nfixed),
                                    as.integer(nvarying),
                                    TRUE ## Remove rows with all missing
                                    ));

        if (numname%in%fixed) {
            while (numname%in%c(fixed)) numname <- paste(numname,"_",sep="")
        }
        if (idname%in%fixed) {
            long <- long[,-(ncol(long)-1)]
            cnames <- c(fixed,vnames,numname)
        } else {
            cnames <- c(fixed,vnames,idname,numname)
        }
        ##  while (idname%in%c(fixed,vnames,numname)) idname <- paste(idname,"_",sep="")
        ##  while (numname%in%c(fixed,vnames)) numname <- paste(numname,"_",sep="")
        
        colnames(long) <- cnames
        if (!numlev) {
            long[,numname] <- factor(long[,numname],labels=thelevels)
        } else {
            if (!identical(order(thelevels),thelevels))
                long[,numname] <- thelevels[long[,numname]]
        }

        if (is_df && factors.keep) { ## Recreate classes 
            vars.orig <- c(fixed,unlist(lapply(varying,function(x) x[1])))
            vars.new <- c(fixed,vnames)
            factors <- which("factor"==classes[vars.orig])
            logicals <- which("logical"==classes[vars.orig])
            for (i in factors) {
                lev <- levels(D0[,vars.orig[i]])
                long[,vars.new[i]] <- factor(lev[long[,vars.new[i]]],levels=lev)
            }
            for (i in logicals) {
                long[,vars.new[i]] <- long[,vars.new[i]]==1
            }
        }
        return(long)
    }

    ## Long to wide format:
    numvar <- idvar <- NULL 
    if (is.character(id)) {
        idvar <- id
        if (length(id)==1) {
            id <- data[,idvar,drop=TRUE]
        } else {
            if (idcombine)
                id <- interaction(as.data.frame(data[,idvar,drop=FALSE]),drop=TRUE)
            else
                id <- data[,idvar[1],drop=TRUE]
        } 
    } else {
        if (length(id)!=nrow(data)) stop("Length of ids and data-set does not agree")
    }
    
    unum <- NULL
    if (!missing(num) && !is.null(num)) {
        if (is.character(num)) {
            numvar <- num
            if (is.character(data[1,num,drop=TRUE])) {
                data[,num,drop=TRUE] <- as.factor(data[,num,drop=TRUE])
            }
            num <- as.integer(data[,num,drop=TRUE])
            if (!labelnum) unum <- sort(unique(data[,numvar,drop=TRUE]))
        } else {
            if (length(num)!=nrow(data)) stop("Length of time and data-set does not agree")
            if (!labelnum) unum <- unique(num)
        }
    } else {
        num <- NULL
    }
    
    nn <- colnames(data)
    if (any(nn=="")) data <- data.frame(data)

    clustud <- cluster.index(id,num=num)
    maxclust <- clustud$maxclust
    idclust <- clustud$idclust  
    obs1 <- clustud$firstclustid+1 ## as.vector(apply(idclust,1,function(x) na.omit(x)[1]))+1

    if (!is.null(numvar)) {
        ii <- which(colnames(data)==numvar)
        data <- data[,-ii,drop=FALSE]
    }

    if (!missing(keep)) {
        keepers <- c(keep,idvar)
        if (!missing(varying)) keepers <- c(keepers,varying)
        ii <- which(colnames(data)%in%keepers)
        data <- data[,ii,drop=FALSE]
    }

    if (missing(varying)) varying <- setdiff(colnames(data),c(idvar))
    vidx <- match(varying,colnames(data))
    N <- nrow(idclust)
    p <- length(varying)
    P <- NCOL(data)  
    fixidx <- setdiff(seq(P),vidx)
    if (is.matrix(data) || (all(apply(data[1,,drop=FALSE],2,is.numeric)) & length(unlist(data[1,]))==length(data[1,]) )) {
        ## Everything numeric - we can work with matrices
        dataw <- matrix(NA, nrow = N, ncol = p * (maxclust-1) + ncol(data))
        dataw[,fixidx] <- as.matrix(data[obs1,fixidx,drop=FALSE])
        mnames <- colnames(data)
        if (!is.null(unum)) {
            mnames[vidx] <- paste(mnames[vidx],unum[1],sep=sep)
        } else {
            mnames[vidx] <- paste(mnames[vidx],1,sep=sep)
        }
        if (p>0) {
            for (i in seq_len(maxclust)) {
                idx <- idclust[, i] + 1
                pos <- vidx
                if (i>1) {
                    pos <- P+seq(p)+p*(i-2)
                }
                dataw[which(!is.na(idx)), pos] <-
                    as.matrix(data[na.omit(idx),vidx,drop=FALSE])      
            }
            if (!is.null(unum)) {
                postn <- unum[-1]
            } else {
                postn <- seq_len(maxclust-1)+1
            }
            ##if (is.null(numname)) postn <- idlev[postn]
            mnames <- c(mnames,
                        as.vector(t(outer(postn,varying,function(x,y) paste(y,x,sep=sep)))))      
        }
        colnames(dataw) <- mnames
        return(dataw)
    }

    ## Potentially slower with data.frame where we use cbind
    for (i in seq_len(maxclust)) {
        if (i==1) {
            dataw <- data[obs1,,drop=FALSE]
            mnames <- names(data);
            dataw[,vidx] <- data[idclust[,i]+1,vidx,drop=FALSE]
            mnames[vidx] <- paste(mnames[vidx],sep,i,sep="")
        } else {
            dataw <- cbind(dataw,data[idclust[,i]+1,varying,drop=FALSE])
            if (!is.null(unum)) 
                mnames <- c(mnames,paste(varying,sep,unum[i],sep=""))
            else
                mnames <- c(mnames,paste(varying,sep,i,sep=""))
        }
    }
    names(dataw) <- mnames    
    return(dataw)
} 



simple.reshape <- function (data, id = "id", num = NULL) {
    cud <- cluster.index(data[, c(id)], num = num, Rindex = 1)
    N <- nrow(cud$idclust)
    p <- ncol(data)
    dataw <- matrix(NA, nrow = N, ncol = p * cud$maxclust)
    for (i in seq_len(cud$maxclust)) {
        dataw[, seq(p) + (i - 1) * p] <- as.matrix(data[cud$idclust[, i] + 1, ])
    }
    colnames(dataw) <- paste(names(data), rep(seq_len(cud$maxclust), each = p), sep = ".")
    return(dataw)
}



