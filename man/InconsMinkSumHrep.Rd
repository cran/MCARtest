% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InconsMinkSumHrep.R
\name{InconsMinkSumHrep}
\alias{InconsMinkSumHrep}
\title{Calculate the H-representation of the general (possibly inconsistent) Minkowski sum}
\usage{
InconsMinkSumHrep(bS, M, round = FALSE)
}
\arguments{
\item{bS}{A binary matrix specifying the set of observation patterns. Each row encodes a single pattern.}

\item{M}{A vector of positive integers giving the alphabet sizes of the discrete variables.}

\item{round}{A logical value indicating whether or not to round coefficients to 15 significant figures.
The function \code{RoundErrors} can be used separately to substitute other values for 15. Defaults to FALSE.}
}
\value{
A halfspace representation object as used by the \code{rcdd} package. See \insertCite{RCDD;textual}{MCARtest} for more detail.
}
\description{
Computes the minimal halfspace representation of the Minkowski sum of the marginal polytope and the inconsistent ball
defined in \insertCite{BS2022;textual}{MCARtest}.
}
\examples{
bS=matrix(c(1,1, 1,0),byrow=TRUE,ncol=2)
InconsMinkSumHrep(bS,c(2,2))
}
\references{
\insertRef{BS2022}{MCARtest}

\insertRef{RCDD}{MCARtest}
}
