#' Computes an inconsistency index for sequences of means.
#'
#' A function that computes the inconsistency index \eqn{M(\mu_\mathbb{S})} for a sequence of
#' means, as defined in Section 5 in \insertCite{BB2024;textual}{MCARtest}.
#'
#' @param mu_S The sequence of means \eqn{\mu_\mathbb{S}}.
#' @param patterns A vector with all the patterns in \eqn{\mathbb{S}}.
#'
#' @references \insertRef{BB2024}{MCARtest}
#'
#' @return The value of \eqn{M()}, in the interval \eqn{[0,1]}.
#' @export
#'
#' @importFrom utils combn
#'
#' @examples
#' library(MASS)
#'
#' d = 3
#' n = 200
#' SigmaS=list() #Random 2x2 correlation matrices (necessarily consistent)
#' for(j in 1:d){
#' x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
#' }
#'
#' X = data.frame(matrix(nrow = 3*n, ncol = 3))
#' X[1:n, c(1,2)] = mvrnorm(n, c(0,0), SigmaS[[1]])
#' X[(n+1):(2*n), c(2, 3)] = mvrnorm(n, c(0,0), SigmaS[[2]])
#' X[(2*n+1):(3*n), c(1, 3)] = mvrnorm(n, c(0,0), SigmaS[[3]])
#' X = as.matrix(X)
#'
#' xxx = get_SigmaS(X)$patterns
#' M(get_SigmaS(X)$muS, xxx)

M = function(mu_S, patterns){
  n_pattern = length(mu_S)
  d = 0
  for (i in 1:n_pattern){
        if (d < max(patterns[[i]])){
            d = max(patterns[[i]])
        }
    }

  max = 0
  if(n_pattern > 1){
    couples = t(combn(1:n_pattern, 2))
    n_couples = dim(t(combn(1:n_pattern, 2)))[1]
    couples = as.list(data.frame(t(combn(1:n_pattern, 2))))

    for (j in 1:d){
      for (i in 1:n_couples){
        if ((j %in% patterns[[ couples$X1[i] ]])&(j %in% patterns[[ couples$X2[i] ]])) {
          candidate = abs(mu_S[[ couples$X1[i] ]][which(patterns[[couples$X1[i]]] == j)[[1]]]
                          - mu_S[[ couples$X2[i] ]][which(patterns[[couples$X2[i]]] == j)[[1]]])
          if (candidate > max){
            max = candidate
          }
        }
      }
    }
  }

  return(max)
}
