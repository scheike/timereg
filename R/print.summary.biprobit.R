##' @export
print.summary.biprobit <- function(x,digits = max(3, getOption("digits") - 2),...) {
  cat("\n")
  printCoefmat(x$coef,digits=digits,...)
  S <- x$score;  names(S) <- rep("",length(S))
  cat("\n")
  print(x$N,quote=FALSE)
  cat("Score: "); cat(formatC(S,...));
  cat("\nlogLik: "); cat(sum(x$logLik),"\n");
  if (!is.null(x$msg)) {
    cat(x$msg,"\n")
  }

  if (!is.null(x$varcomp)) {
    cat("\n")
    res <- x$varcomp
    if (!is.null(x$prob)) {
      res <- rbind(res,x$prob)
    }
    idx <- unlist(sapply(c("Concordance","Marginal"),function(x) grep(x,rownames(res))))
    idx2 <- setdiff(seq(nrow(res)),idx)
    res2 <- rbind(res[idx2,],rep(NA,ncol(res)),res[idx,])
    print(RoundMat(res2,digits=digits,na=FALSE),quote=FALSE)
  }
  cat("\n")
}
