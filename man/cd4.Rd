\name{cd4}
\alias{cd4}
\non_function{}
\title{The multicenter AIDS cohort study
}
\description{CD4 counts collected over time.}
\format{
  This data frame contains the following columns:
\describe{
    \item{obs}{a numeric vector. Number of observations.}
    \item{id}{a numeric vector. Id of subject.}
    \item{visit}{
      a numeric vector. Timings of the visits in years.}
    \item{smoke}{a numeric vector code. 0: non-smoker, 1: smoker.}
    \item{age}{a numeric vector. Age of the patient at the start of the trial.}
    \item{cd4}{a numeric vector. CD4 percentage at the current visit.}
    \item{cd4.prev}{a numeric vector. CD4 level at the preceding visit.}
    \item{precd4}{a numeric vector. Post-infection CD4 percentage.}
    \item{lt}{a numeric vector. Gives the starting time for the
      time-intervals.}
    \item{rt}{a numeric vector. Gives the stopping time for the
      time-interval.}
}
}
\source{
 MACS Public Use Data Set Release PO4 (1984-1991). See reference.
}
\references{
Kaslow et al. (1987), The multicenter AIDS cohort study: rational,
organisation and selected characteristics of the participants. 
Am. J. Epidemiology 126, 310--318. 
}
\examples{
data(cd4)
names(cd4)
}
\keyword{datasets}
