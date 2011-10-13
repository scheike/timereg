.packageName <- "compRisk"

.First.lib <- function(lib, pkg) {
  library.dynam("compRisk", pkg, lib)
}

.Last.lib <- function(lib){
  library.dynam.unload("compRisk",lib)
}

