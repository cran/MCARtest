% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Amatrix.R
\name{Amatrix}
\alias{Amatrix}
\title{Generate the matrix A, whose columns are the vertices of the marginal polytope.}
\usage{
Amatrix(bS, M)
}
\arguments{
\item{bS}{A binary matrix specifying the set of observation patterns. Each row encodes a single pattern.}

\item{M}{A vector of positive integers giving the alphabet sizes of the discrete variables.}
}
\value{
The matrix A.
}
\description{
Generate the matrix A, whose columns are the vertices of the marginal polytope.
}
\examples{
bS=matrix(c(1,1,0, 1,0,1, 0,1,1),byrow=TRUE,ncol=3)
M=c(2,2,2)
Amatrix(bS,M)
}
