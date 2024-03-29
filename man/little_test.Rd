% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/little_test.R
\name{little_test}
\alias{little_test}
\title{Carry out Little's test of MCAR}
\usage{
little_test(X, alpha, type = "mean&cov")
}
\arguments{
\item{X}{The dataset with incomplete data, where all the pairs of variables are observed together.}

\item{alpha}{The nominal level of the test.}

\item{type}{Determines the test statistic to use, based on the discussion in Section 5 in \insertCite{BB2024;textual}{MCARtest}.
The default option is "mean&cov", and uses the test statistic \eqn{d^2_{\mathrm{aug}}}. When set equal to "cov", implements a test
of MCAR based on \eqn{d^2_{\mathrm{cov}}}, while, when set equal to "mean", implements the classical Little's test as defined in
\insertCite{Little1988;textual}{MCARtest}.}
}
\value{
A Boolean, where TRUE stands for reject MCAR. This is computed by comparing the p-value of Little's test,
found by comparing the log likelihood ratio statistic to the chi-squared distribution with the appropriate number
of degrees of freedom, with the nominal level \code{alpha}. Described in \insertCite{Little1988;textual}{MCARtest}.
}
\description{
Carry out Little's test of MCAR
}
\examples{
library(MASS)
alpha = 0.05
n = 200

SigmaS=list() #Random 2x2 correlation matrices (necessarily consistent)
for(j in 1:3){
x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1)
SigmaS[[j]]=cov2cor(x\%*\%t(x) + y\%*\%t(y))
}

X1 = mvrnorm(n, c(0,0), SigmaS[[1]])
X2 = mvrnorm(n, c(0,0), SigmaS[[2]])
X3 = mvrnorm(n, c(0,0), SigmaS[[3]])
columns = c("X1","X2","X3")
X = data.frame(matrix(nrow = 3*n, ncol = 3))
X[1:n, c("X1", "X2")] = X1
X[(n+1):(2*n), c("X2", "X3")] = X2
X[(2*n+1):(3*n), c("X1", "X3")] = X3
X = as.matrix(X)

little_test(X, alpha)
}
\references{
\insertRef{BB2024}{MCARtest}

\insertRef{Little1988}{MCARtest}
}
