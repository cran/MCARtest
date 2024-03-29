% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Csimple.R
\name{Csimple}
\alias{Csimple}
\title{Calculate the critical value for our simple test}
\usage{
Csimple(nS, bS, M, alpha)
}
\arguments{
\item{nS}{A vector of sample sizes, with each entry corresponding to an observation pattern.}

\item{bS}{A binary matrix specifying the set of observation patterns. Each row encodes a single pattern.}

\item{M}{A vector of positive integers giving the alphabet sizes of the discrete variables.}

\item{alpha}{The desired significance level \eqn{\alpha} of the test.}
}
\value{
The universal critical value defined in \insertCite{BS2022;textual}{MCARtest}.
}
\description{
Calculate a simple critical value for an MCAR test using only knowledge of
the set of observation patterns and the joint observation space.
}
\examples{
bS=matrix(c(1,1,0, 1,0,1, 0,1,1),byrow=TRUE,ncol=3)
r=4; s=3
M=c(r,s,2)
Csimple(rep(1000,3),bS,M,0.05)
}
\references{
\insertRef{BS2022}{MCARtest}
}
