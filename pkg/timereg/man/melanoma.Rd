\name{melanoma}
\alias{melanoma}
\non_function{}
\title{The Melanoma Survival Data}
\description{
  The melanoma data frame has 205 rows and 7 columns.
  It contains data relating to survival of patients after 
  operation for malignant melanoma collected at 
  Odense University Hospital by K.T.  Drzewiecki. 
}
\format{
This data frame contains the following columns:
\describe{
\item{no}{
a numeric vector. Patient code.
}
\item{status}{
a numeric vector code. Survival status. 1: dead from melanoma, 2: alive,
3: dead from other cause.
}
\item{days}{
a numeric vector. Survival time.
}
\item{ulc}{
a numeric vector code. Ulceration, 1: present, 0: absent.
}
\item{thick}{
a numeric vector. Tumour thickness (1/100 mm).
}
\item{sex}{
a numeric vector code. 0: female, 1: male.
}
}
}
\source{
  Andersen, P.K., Borgan O, Gill R.D., Keiding N. (1993),
  \emph{Statistical Models Based on 
  Counting Processes}, Springer-Verlag. 

  Drzewiecki, K.T., Ladefoged, C., and Christensen, H.E. (1980),
  Biopsy and prognosis for cutaneous malignant melanoma in clinical
  stage I. Scand. J. Plast. Reconstru. Surg. 14, 141-144. 
}
\examples{
data(melanoma)
names(melanoma)
}
\keyword{datasets}
