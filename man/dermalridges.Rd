\docType{data}
\name{dermalridges}
\alias{dermalridges}
\title{Dermal ridges data (families)}
\format{Data on 50 families with ridge counts in left and right
hand for moter, father and each child. Family id in 'family' and
gender and child number in 'sex' and 'child'.}
\source{
Sarah B. Holt (1952). Genetics of dermal ridges: bilateral
asymmetry in finger ridge-counts.  Annals of Eugenics 17 (1),
pp.211--231. DOI: 10.1111/j.1469-1809.1952.tb02513.x
}
\description{
Data on dermal ridge counts in left and right hand in (nuclear) families
}
\examples{
data(dermalridges)
fast.reshape(dermalridges,id="family",varying=c("child.left","child.right","sex"))
}
\keyword{data}

