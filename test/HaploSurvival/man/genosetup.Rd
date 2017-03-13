\name{geno.setup}
\alias{geno.setup}
\alias{print.geno.setup}
\alias{summary.geno.setup}
\title{Infers possible haplotypes based on genotype data.}
\description{

Input is genotype data. 

Can handle all genotype data with multialleles and is thus not restricted to
SNP data. 
}
\usage{
geno.setup(G,returnNames=TRUE,haplo.baseline=NULL)
}
\arguments{
\item{G}{genotype data, n rows with 2*p columns. }
\item{returnNames}{return names of haplotype.}
\item{haplo.baseline}{gives the name of the haplotype that should be used as 
baseline, for latter MLE fitting.  Default is last haplotype.}
}
\details{
  This function provides an efficient representation of haplotypes that
  are consistent with genotype data. This speeds up the calculation of
  the likelihood of the genotype data, given the frequencies of
  haplotypes.

  The component \code{out$HPIordered} returned contains the most important
  information. \code{out$HPIordered} is list, with one entry per
  person. \code{out$HPIordered[[i]]} is a list, with one entry per
  haplotype pair consistent with the genotype data
  \code{out$HPIordered[[i]][[j]]} is a vector of two integers, which are
  the indices of the \code{j}th pair of haplotypes consistent with the
  genotypes for
  individual \code{i}. The haplotypes that these indices correspond to
  can be examined by looking at 

  \code{out$uniqueHaploNames[out$HPIordered[[i]][[j]]]}.
}
\value{
Returns an object of class \code{geno.setup} with the following arguments:
\item{HPIordered}{list of indices of haplotype pairs that can lead to the 
observed genotype}
\item{uniqueHaploNames}{names of haplotypes, with alleles separated by
  columns. This list includes only haplotypes which are consistent with
  the genotypes for individuals in the sample.}
\item{nAllelesPerLocus}{number of alleles for each locus.} 
\item{unorderedAlleles}{A list, with one component corresponding to a
  locus, and this component is a vector of the allele names found at
  this locus.}
\item{nPeople}{number of subjects.} 
\item{nPossHaps}{A vector of the number of possible haplotype pairs for each subject.}
\item{nLoci}{number of loci.}
\item{initial.freqs}{
  simple estimator of haplotype frequencies, based on simply counting
  the occurrences of each haplotype relative to the total number of
  possible haplotypes for all subjects.}
}
\author{Jeremy Silver}
\examples{
n<-10
g<-matrix(rbinom(10*4,1,0.5),n,4); # simple genotypes

gs<-geno.setup(g); 
summary(gs); 

gs
}
\keyword{survival}
