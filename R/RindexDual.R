#' A function computing the incompatibility index and associated closest joint mass function using the dual formulation
#'
#' A function solving a linear program to compute the incompatibility index \eqn{R()} defined in \insertCite{BS2022;textual}{MCARtest},
#' in the case of having discrete random variables.
#' Uses \code{Amatrix} to define to constraint matrix and \code{lpSolve} to implement the linear optimisation.
#'
#' @param pS A sequence of probability mass functions on the marginal spaces.
#' @param bS A binary matrix specifying the set of observation patterns. Each row encodes a single pattern.
#' @param M A vector of positive integers giving the alphabet sizes of the discrete variables.
#' @param lp_solver An argument passed to HiGHS specifying which solver to use. See \insertCite{highs;textual}{MCARtest} for more detail.
#' @param simplex_strategy An argument passed to HiGHS specifying which solver to use. See \insertCite{highs;textual}{MCARtest} for more detail.
#'
#' @return The value of \eqn{R()}, in the interval \eqn{[0,1]}.
#' @return The optimal solution to the linear program
#' @export
#'
#' @references
#' \insertRef{BS2022}{MCARtest}
#' \insertRef{highs}{MCARtest}
#'
#' @importFrom highs highs_solve
#'
#' @examples
#' bS=matrix(c(1,1,0, 1,0,1, 0,1,1),byrow=TRUE,ncol=3)
#' M=c(2,2,2)
#' A=Amatrix(bS,M)
#'
#' pS=rep(0.25,12)
#' linprog=RindexDual(pS,bS,M)
#' rbind(pS,as.vector(A%*%linprog[[2]])/(1-linprog[[1]]))
#'
#' pS=c(0.125,0.375,0.375,0.125,0.250,0.250,0.250,0.250,0.100,0.400,0.400,0.100)
#' linprog=RindexDual(pS,bS,M)
#' rbind(pS,as.vector(A%*%linprog[[2]])/(1-linprog[[1]]))
#'

RindexDual <- function(pS, bS, M, lp_solver = "default", simplex_strategy = 4) {
  totCardS <- length(pS)
  cardChi <- prod(M)

  A <- AmatrixSparse(bS, M)

  if (lp_solver != "default") {
    Sol <- highs::highs_solve(L = rep(1, cardChi), lower = rep(0, cardChi),
                              upper = rep(Inf, cardChi), A = A, lhs = rep(-Inf, totCardS),
                              rhs = pS, types = rep(1, cardChi), maximum = TRUE,
                              control = highs::highs_control(log_to_console = FALSE, solver = lp_solver, simplex_strategy = simplex_strategy))

    if (Sol$status == 4) {
      return(list("Failed", "Failed"))
    }

    return(list(1 - Sol$objective_value, Sol$primal_solution))
  }

  Sol <- highs::highs_solve(L = rep(1, cardChi), lower = rep(0, cardChi),
                            upper = rep(Inf, cardChi), A = A, lhs = rep(-Inf, totCardS),
                            rhs = pS, types = rep(1, cardChi), maximum = TRUE,
                            control = highs::highs_control(log_to_console = FALSE, solver = "ipm"))
  if (Sol$status == 4) {
    Sol <- highs::highs_solve(L = rep(1, cardChi), lower = rep(0, cardChi),
                              upper = rep(Inf, cardChi), A = A, lhs = rep(-Inf, totCardS),
                              rhs = pS, types = rep(1, cardChi), maximum = TRUE,
                              control = highs::highs_control(log_to_console = FALSE, solver = "simplex", simplex_strategy = simplex_strategy))
  }

  return(list(1 - Sol$objective_value, Sol$primal_solution))
}
