##' @export
sim.cox <- function(x,...) {
    timereg::sim.cox(cox=x,...)
}

##' @export
simulate.cox <- function(object,...) sim(object,...)
