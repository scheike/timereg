\name{mela.pop}
\alias{mela.pop}
\non_function{}
\title{Melanoma data and Danish population mortality by age and sex
}
\description{Melanoma data with background mortality of Danish population.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{id}{ 
      a numeric vector. Gives patient id.
    }
    \item{sex}{
      a numeric vector. Gives sex of patient.
    }
    \item{start}{
      a numeric vector. 
Gives the starting time for the
      time-interval for which the covariate rate is representative.
    }
    \item{stop}{
      a numeric vector. Gives the stopping time for the
      time-interval for which the covariate rate is representative.
    }
\item{status}{
a numeric vector code. Survival status. 1: dead from melanoma, 0: alive or
 dead from other cause.
}
    \item{age}{
      a numeric vector. Gives the age of the patient at removal of tumor.
    }
    \item{rate}{
      a numeric vector. Gives the population mortality for the 
                   given sex and age.  Based on Table A.2 in Andersen 
                   et al. (1993). 
    }
  }
}
\source{
  Andersen, P.K., Borgan O, Gill R.D., Keiding N. (1993),
  \emph{Statistical Models Based on 
  Counting Processes}, Springer-Verlag. 
}
\examples{
data(mela.pop)
names(mela.pop)
}
\keyword{datasets}
