#' Compute the columnwise average of means/variances
#'
#' A function that computes \eqn{\bar{\operatorname{av}}_j(\mu_{\mathbb{S}})} as defined in
#' Section 5 in \insertCite{BB2024;textual}{MCARtest}, or
#' \eqn{\bar{\operatorname{av}}_j(\sigma^2_{\mathbb{S}})} as defined in Section 2 in
#' \insertCite{BB2024;textual}{MCARtest}. The sequence of means/variances, and the
#' sequence of patterns, are calculated with \code{getSigmaS}.
#'
#' @param type If set equal to "mean", computes \eqn{\bar{\operatorname{av}}_j(\mu_{\mathbb{S}})}.
#' If set equal to "var", computes \eqn{\bar{\operatorname{av}}_j(\sigma^2_{\mathbb{S}})}.
#' @param X The whole dataset with missing values.
#'
#' @references \insertRef{BB2024}{MCARtest}
#'
#' @return The value of \eqn{\bar{\operatorname{av}}_j(\sigma^2_{\mathbb{S}})} or
#' \eqn{\bar{\operatorname{av}}_j(\mu_{\mathbb{S}})}.
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
#' compute_av("var", X)
#' compute_av("mean", X)

compute_av = function(type, X){

  result = get_SigmaS(X)
  d = result$ambient_dimension

  if (type=="mean"){
    value_S = result$muS
  }
  else if (type=="var"){
    value_S = result$sigma_squared_S
  }
  else {
    print("An error occured")
  }

  patterns = result$pattern
  card_patterns = length(patterns)
  n_pattern = result$n_pattern
  data_pattern = result$data_pattern

  res = numeric(length = d)
  for (j in 1:d){
    tmp = 0
    card_Sj = 0
    for (i in 1:n_pattern){
      if (j %in% patterns[[i]]){
        tmp = tmp + value_S[[i]][which(patterns[[i]] == j)[[1]]]
        card_Sj = card_Sj+1
      }
    }

    res[j] = tmp/card_Sj
  }

  return(res)
}
