% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EquivalenceClass.R
\name{EquivalenceClass}
\alias{EquivalenceClass}
\title{Simplifies H-representation by exploiting symmetry}
\usage{
EquivalenceClass(bS, M, Hrep)
}
\arguments{
\item{bS}{A binary matrix specifying the set of observation patterns. Each row encodes a single pattern.}

\item{M}{A vector of positive integers giving the alphabet sizes of the discrete variables.}

\item{Hrep}{An H-representation generated by \code{MargPolyHrep}, \code{ConsMinkSumHrep} or \code{InconsMinkSumHrep}.}
}
\value{
A list of representative halfspace normals.
}
\description{
The marginal polytope and related objects have many symmetries. By relabelling the levels
of discrete variables we transform facets into other facets. This function reduces a list of
halfspace normals to its equivalence classes.
}
\examples{
bS=matrix(c(1,1,0, 1,0,1, 0,1,1),byrow=TRUE,ncol=3) # Our canonical 3d example
Hrep=MargPolyHrep(bS,c(2,2,2))
EquivalenceClass(bS,c(2,2,2),Hrep)

}
