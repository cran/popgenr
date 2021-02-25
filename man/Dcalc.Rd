\name{Dcalc}
\alias{Dcalc}
\title{Calculate and visualize LD}
\description{
  Calculates the degree of linkage disequilibirum between two biallelic diploid loci and plots the results.
}
\usage{
Dcalc(AABB, AaBB, aaBB, AABb, AaBb, aaBb, AAbb, Aabb, aabb)
}
\arguments{
  \item{AABB}{Count of A and B double homozygotes.}
  \item{AaBB}{Count of A heterozygotes, B homozygotes.}
  \item{aaBB}{Count of a and B double homozygotes.}
  \item{AABb}{Count of A homozygotes, B heterozygotes.}
  \item{AaBb}{Count of double heterozygotes.}
  \item{aaBb}{Count of a homozygotes, B heterozygote.}
  \item{AAbb}{Count of A and b double homozygotes.}
  \item{Aabb}{Count of A heterozygotes, b homozygotes.}
  \item{aabb}{Count of a and b double homozygotes.}
}

\examples{
#assign genotype counts
AABB=2
AaBB=0
aaBB=0
AABb=0
AaBb=1
aaBb=0
AAbb=1
Aabb=0
aabb=0

#run function
Dcalc(AABB, AaBB, aaBB, AABb, AaBb, aaBb, AAbb, Aabb, aabb)
}
\keyword{file}
