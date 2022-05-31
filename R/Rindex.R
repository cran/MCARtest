#' A function computing the incompatibility index
#'
#' A function solving a linear program to compute the incompatibility index \eqn{R()} defined in \insertCite{BS2022;textual}{MCARtest},
#' in the case of having discrete random variables.
#' Uses \code{Amatrix} to define to constraint matrix and \code{lpSolve} to implement the linear optimisation.
#'
#' @param pS A sequence of probability mass functions on the marginal spaces.
#' @param bS A binary matrix specifying the set of observation patterns. Each row encodes a single pattern.
#' @param M A vector of positive integers giving the alphabet sizes of the discrete variables.
#'
#' @return The value of \eqn{R()}, in the interval \eqn{[0,1]}.
#' @export
#'
#' @references
#' \insertRef{BS2022}{MCARtest}
#'
#' @importFrom lpSolve lp
#'
#' @examples
#' bS=matrix(c(1,1,0, 1,0,1, 0,1,1),byrow=TRUE,ncol=3)
#' M=c(2,2,2)
#'
#' pS=rep(0.25,12)
#' Rindex(pS,bS,M)
#'
#' pS=c(0.125,0.375,0.375,0.125,0.250,0.250,0.250,0.250,0.100,0.400,0.400,0.100)
#' Rindex(pS,bS,M)

Rindex=function(pS,bS,M){
  A=Amatrix(bS,M)
  linprog=lp(direction="min",
             objective.in = pS,
             const.mat = A,
             transpose.constraints = FALSE,
             const.dir = rep(">=",prod(M)),
             const.rhs = rep(1,prod(M))
  )
  return(1-linprog$objval)
}
