print.CRR <- function(x,...){
  f <- x$crrFit
  f$call <- x$call
  print(summary(f,...))
  cat(paste("\nConvergence:",f$converged),"\n\n")
}
