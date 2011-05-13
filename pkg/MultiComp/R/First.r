.packageName <- "MultiComp"

.First.lib <- function(lib, pkg) {
  library.dynam("MultiComp", pkg, lib)
  cat("This is MultiComp \n\n");
}

.Last.lib <- function(lib){
  library.dynam.unload("MultiComp",lib)
}
