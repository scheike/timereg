##' @S3method print summary.biprobit
print.summary.biprobit <- function(x,digits = max(3, getOption("digits") - 2),...) {
  cat("\n")
  printCoefmat(x$coef,digits=digits,...)
##  print(x$object,digits=digits)
  S <- colSums(x$score);  names(S) <- rep("",length(S))
  cat("\n")
  print(x$N,quote=FALSE)
  ##  suppressMessages(browser())
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
    print(RoundMat(res,digits=digits),quote=FALSE)
  }
  cat("\n")
}
