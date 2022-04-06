\name{TRACE}
\alias{TRACE}
\alias{sTRACE}
\alias{tTRACE}
\non_function{}
\title{The TRACE study group of myocardial infarction}
\description{
  The TRACE data frame contains 1877 patients and is a subset of
  a data set consisting of approximately 6000 patients. 
  It contains data relating survival of patients after 
  myocardial infarction to various risk factors.

  sTRACE is a subsample consisting of 300 patients.

  tTRACE is a subsample consisting of 1000 patients.
}
\format{
This data frame contains the following columns:
\describe{
\item{id}{a numeric vector. Patient code.
}
\item{status}{
a numeric vector code. Survival status. 9: dead from myocardial infarction, 0: alive, 7,8: dead from other causes. 
}
\item{time}{
a numeric vector. Survival time in years.
}
\item{chf}{
a numeric vector code. Clinical heart pump failure, 1: present, 0: absent.
}
\item{diabetes}{
a numeric vector code. Diabetes, 1: present, 0: absent.
}
\item{vf}{
a numeric vector code. Ventricular fibrillation, 1: present, 0: absent.
}
\item{wmi}{
a numeric vector.  Measure of heart pumping effect based on ultrasound
measurements where 2 is normal and 0 is worst.
}
\item{sex}{
a numeric vector code. 1: female, 0: male.
}
\item{age}{
a numeric vector code. Age of patient.
}
}
}
\source{
The TRACE study group. 

Jensen, G.V., Torp-Pedersen, C., Hildebrandt, P., Kober, L.,
F. E. Nielsen, Melchior, T., Joen, T. and P. K. Andersen (1997),
Does in-hospital ventricular fibrillation
affect prognosis after myocardial  infarction?,
European Heart Journal 18, 919--924.
}
\examples{
data(TRACE)
names(TRACE)
}
\keyword{datasets}
