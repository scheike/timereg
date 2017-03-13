\name{diabetes}
\alias{diabetes}
\non_function{}
\title{The Diabetic Retinopathy Data}
\description{
The data was colleceted to test a laser treatment for delaying
blindness in patients with dibetic retinopathy. The subset of 197
patiens given in Huster et al. (1989) is used. 
}
\format{
This data frame contains the following columns:
\describe{
\item{id}{a numeric vector. Patient code.}
\item{agedx}{a numeric vector. Age of patient at diagnosis.}
\item{time}{a numeric vector. Survival time: time to blindness or 
censoring.}
\item{status}{
a numeric vector code. Survival status. 1: blindness, 0: censored.}
\item{trteye}{a numeric vector code. Random eye selected for
treatment. 1: left eye 2: right eye.}
\item{treat}{a numeric vector. 1: treatment 0: untreated.} 
\item{adult}{a numeric vector code. 1: younger than 20, 2: older than 20.}
}
}
\source{
Huster W.J. and Brookmeyer, R. and Self. S. (1989) 
MOdelling paired survival data with covariates, Biometrics 45,
145-56.
}
\examples{
data(diabetes)
names(diabetes)
}
\keyword{datasets}
