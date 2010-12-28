.packageName <- "timereg"

.First.lib <- function(lib, pkg) {
library.dynam("timereg", pkg, lib)
}

.Last.lib <- function(lib){
  library.dynam.unload("timereg",lib)
}
