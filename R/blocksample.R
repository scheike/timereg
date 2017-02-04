##' @export
dsample <- function(data,x,size=NULL,replace=TRUE,...) {
    if (missing(x)) {
        if (is.null(size)) size <- NROW(data)
        return(data[sample.int(NROW(data), size, replace=replace),,drop=FALSE])
    }
    inp <- procform(x,data=data,return.formula=TRUE)
    if (length(inp$filter.expression)>0)
        data <- subset(data,eval(inp$filter.expression))
    if (!is.null(inp$predictor))
        idvar <- model.frame(inp$predictor, data=data, na.action=na.pass,...)
    if (!is.null(inp$response)) {
        if (is.null(size)) size <- NROW(data)
        data <- model.frame(inp$response, data=data, na.action=na.pass,...)
    }
    if (is.null(inp$predictor)) {
        if (is.null(size)) size <- NROW(data)
        return(data[sample.int(NROW(data), size, replace=replace),,drop=FALSE])
    }
    blocksample(data,idvar=idvar,size=size,replace=replace,...)
}

##' Sample blockwise from clustered data
##'
##' @title Block sampling
##' @param data Data frame
##' @param idvar Column defining the clusters
##' @param size Size of samples
##' @param replace Logical indicating wether to sample with replacement
##' @param \dots additional arguments to lower level functions
##' @return \code{data.frame}
##' @author Klaus K. Holst
##' @keywords models utilities
##' @aliases blocksample dsample
##' @details Original id is stored in the attribute 'id'
##' @export
##' @examples
##'
##' d <- data.frame(x=rnorm(5), z=rnorm(5), id=c(4,10,10,5,5), v=rnorm(5))
##' (dd <- blocksample(d,size=20,~id))
##' attributes(dd)$id
##'
##' \dontrun{
##' blocksample(data.table::data.table(d),1e6,~id)
##' }
##'
##'
##' d <- data.frame(x=c(1,rnorm(9)),
##'                z=rnorm(10),
##'                id=c(4,10,10,5,5,4,4,5,10,5),
##'                id2=c(1,1,2,1,2,1,1,1,1,2),
##'                v=rnorm(10))
##' dsample(d,~id, size=2)
##' dsample(d,.~id+id2)
##' dsample(d,x+z~id|x>0,size=5)
##'
blocksample <- function(data, size, idvar=NULL, replace=TRUE, ...) {
    if (is.null(idvar)) {
        return(data[sample(NROW(data),size,replace=replace),,drop=FALSE])
    }
    if (inherits(idvar,"formula")) {
        idvar <- all.vars(idvar)
    }
    if (NROW(idvar)==nrow(data)) {
        id0 <- idvar
    } else {
        if (inherits(data,"data.table")) {
            id0 <- as.data.frame(data[,idvar,with=FALSE])[,1]
        } else id0 <- data[,idvar]
    }
    if (NCOL(id0)>1) {
        id0 <- interaction(as.data.frame(id0))
    }
    id0 <- as.matrix(id0)[,1,drop=TRUE]
    ii <-  cluster.index(as.matrix(id0))
    size <- ifelse(missing(size) || is.null(size),ii$uniqueclust,size)
    ids <- sample(seq(ii$uniqueclust), size=size,replace=replace)
    idx <- na.omit(as.vector(t(ii$idclustmat[ids,])))+1
    newid <- rep(seq(size), ii$cluster.size[ids])
    oldid <- id0[idx]
    res <- data[idx,]
    if (is.character(idvar) && length(idvar)==1) {
        res[,idvar] <- newid
    } else {
        res <- cbind(res,id=newid)
        colnames(res) <- make.unique(colnames(res))
    }
    attributes(res)$id <- oldid
    return(res)
}
