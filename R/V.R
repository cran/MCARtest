#' Computes an inconsistency index for sequences of variances.
#'
#' A function that computes the inconsistency index \eqn{V(\sigma^2_\mathbb{S})} for a sequence of
#' variances as defined in Section 2 in \insertCite{BB2024;textual}{MCARtest}, given the fact that
#' \eqn{\bar{\operatorname{av}}(\sigma^2_{\mathbb{S}_j}) = 1}.
#'
#' @param sigma_squared_S The sequence of variances \eqn{\sigma_\mathbb{S}^2}.
#' @param patterns A vector with all the patterns in \eqn{\mathbb{S}}.
#'
#' @references \insertRef{BB2024}{MCARtest}
#'
#' @return The value of \eqn{V()}, in the interval \eqn{[0,1]}.
#' @export
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
#' av_sigma = compute_av("var", X)
#' X_new = X
#' for (j in 1:3){
#'   X_new[,j] = X[,j]/sqrt(av_sigma[j])
#'   }
#'
#' V(get_SigmaS(X_new)$sigma_squared_S, xxx)

V = function(sigma_squared_S, patterns){
  n_pattern = length(sigma_squared_S)
  d = 0
  for (i in 1:n_pattern){
    if (d < max(patterns[[i]])){
      d = max(patterns[[i]])
    }
  }

  min = 1

  for (j in 1:d){
    for (i in 1:n_pattern){
      if (j %in% patterns[[i]]) {
        candidate = sigma_squared_S[[i]][which(patterns[[i]] == j)[[1]]]
        if (candidate < min){
          min = candidate
        }
      }
    }
  }

  return(1-min)
}
