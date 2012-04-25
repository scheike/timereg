##' @S3method print summary.biprobit
print.summary.biprobit <- function(x,digits = max(3, getOption("digits") - 2),...) {
  cat("\n")
  print(x$object,digits=digits)
  S <- colSums(x$object$score);  names(S) <- rep("",length(S))
  cat("\n")
  print(x$object$N,quote=FALSE)
  ##  suppressMessages(browser())
  cat("Score: "); cat(formatC(S,...));
  cat("\nlogLik: "); cat(sum(x$object$logLik),"\n");
  if (!is.null(x$object$msg)) {
    cat(x$object$msg,"\n")
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
