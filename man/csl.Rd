\name{csl}
\alias{csl}
\non_function{}
\title{CSL liver chirrosis data}
\description{Survival status for the liver chirrosis patients of Schlichting et al.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{id}{
      a numeric vector. Id of subject.
    }
    
    \item{time}{
      a numeric vector. Time of measurement.
    }
 \item{prot}{
      a numeric vector. Prothrombin level at measurement time.
    }
   \item{dc}{
      a numeric vector code. 0: censored observation, 1: died at
      eventT.
    }
    \item{eventT}{
      a numeric vector. Time of event (death).
    }
    \item{treat}{
      a numeric vector code. 0:  active treatment
      of prednisone, 1: placebo treatment.
    }
    \item{sex}{
      a numeric vector code. 0: female, 1: male.
    }
    \item{age}{
      a numeric vector. Age of subject at inclusion time subtracted 60.
    }
      \item{prot.base}{
      a numeric vector. Prothrombin base level before
      entering the study.
    }
    \item{prot.prev}{
      a numeric vector. Level of prothrombin at previous measurement
      time.
    }
    \item{lt}{
      a numeric vector. Gives the starting time for the
      time-intervals.
    }
    \item{rt}{
      a numeric vector. Gives the stopping time for the
      time-intervals.
    }
      }
}
\source{P.K. Andersen}
\references{
Schlichting, P., Christensen, E., Andersen, P., Fauerholds, L.,
Juhl, E., Poulsen, H. and  Tygstrup, N. (1983), The Copenhagen 
Study Group for Liver Diseases, Hepatology 3, 889--895
}
\examples{
data(csl)
names(csl)
}
\keyword{datasets}
