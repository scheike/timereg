\name{summary.haplo.freqs}
\alias{summary.haplo.freqs}
\title{Prints summary statistics for haplo frequency fit}
\description{
Presents MLE of haplotype frequencies based on Hardy-Weinberg 
equilibrium.
}
\usage{
\method{summary}{haplo.freqs}(object,...)
}
\arguments{
\item{object}{an haplo.freq object.}
\item{...}{unused arguments - for S3 compatibility}
}
\references{
}
\author{Jeremy Silver}
\examples{
### Simulating data
n<-200;
haplo.types<-rbind(c(0,0),c(1,0),c(0,1),c(1,1))
pairs<-matrix(sample(1:4,n*2,replace=TRUE),n,2) ### randomly selects haplotypes

g<-matrix(haplo.types[pairs,],n,4) ### observed genotype, haplotype-pair 1 = g[c(1,3),]

gs<-geno.setup(g);

hapfit<-haplo.freqs(g,geno.setup=gs)
summary(hapfit); 

X<-matrix(1,3,1);
hapfit<-haplo.freqs(g,geno.setup=gs,haplo.design=X)
summary(hapfit); 
}
\keyword{survival}
