\name{coal}
\alias{coal}
\title{Simulate and visualize a coalescent process}
\description{
  A function that provides a simple starting off point to simulate a coalescent process.
}
\usage{
coal(length, number, muscale, reps, prnt)
}
\arguments{
  \item{length}{Length of sequence to simulate.}
  \item{number}{Number of starting lineages.}
  \item{muscale}{Mutation rate, scaled by 4Neu.}
  \item{reps}{Number of replicates to run.}
  \item{prnt}{Whether to print calculation output to the terminal (prnt=1), or not (prnt=0).}
}

\examples{
#assign genotype counts
length=250

number=8

muscale=6.25 # theta = 4Neu = 25, theta/4 = 6.25

reps=100

prnt=0

#run function
coal(length, number, muscale, reps, prnt)
}
\keyword{file}
