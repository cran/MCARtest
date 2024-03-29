% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RindexDual.R
\name{RindexDual}
\alias{RindexDual}
\title{A function computing the incompatibility index and associated closest joint mass function using the dual formulation}
\usage{
RindexDual(pS, bS, M, lp_solver = "default", simplex_strategy = 4)
}
\arguments{
\item{pS}{A sequence of probability mass functions on the marginal spaces.}

\item{bS}{A binary matrix specifying the set of observation patterns. Each row encodes a single pattern.}

\item{M}{A vector of positive integers giving the alphabet sizes of the discrete variables.}

\item{lp_solver}{An argument passed to HiGHS specifying which solver to use.}

\item{simplex_strategy}{An argument passed to HiGHS specifying which solver to use.}
}
\value{
The value of \eqn{R()}, in the interval \eqn{[0,1]}.

The optimal solution to the linear program
}
\description{
A function solving a linear program to compute the incompatibility index \eqn{R()} defined in \insertCite{BS2022;textual}{MCARtest},
in the case of having discrete random variables.
Uses \code{Amatrix} to define to constraint matrix and \code{lpSolve} to implement the linear optimisation.
}
\examples{
bS=matrix(c(1,1,0, 1,0,1, 0,1,1),byrow=TRUE,ncol=3)
M=c(2,2,2)
A=Amatrix(bS,M)

pS=rep(0.25,12)
linprog=RindexDual(pS,bS,M)
rbind(pS,as.vector(A\%*\%linprog[[2]])/(1-linprog[[1]]))

pS=c(0.125,0.375,0.375,0.125,0.250,0.250,0.250,0.250,0.100,0.400,0.400,0.100)
linprog=RindexDual(pS,bS,M)
rbind(pS,as.vector(A\%*\%linprog[[2]])/(1-linprog[[1]]))

}
\references{
\insertRef{BS2022}{MCARtest}
}
