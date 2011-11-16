.packageName <- "MultiComp"

.First.lib <- function(lib, pkg) {
  library.dynam("MultiComp", pkg, lib)
}

.Last.lib <- function(lib){
  library.dynam.unload("MultiComp",lib)
}
