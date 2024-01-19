#' Computes the sequence of patterns, means, variances, covariance and correlation
#' matrices for a given dataset with missing values.
#'
#' Using the same the notation of  \insertCite{BB2024;textual}{MCARtest}, computes
#' the sequence of patterns \eqn{\mathbb{S}}, means \eqn{\mu_\mathbb{S}}, variances
#' \eqn{\sigma^2_\mathbb{S}}, and correlation matrices \eqn{\Sigma_\mathbb{S}}
#' for a dataset with missing values.
#'
#' @param X The dataset with incomplete data.
#'
#'
#' @returns \code{patterns} The sequence of patterns \eqn{\mathbb{S}}.
#' @returns \code{n_pattern} The cardinality of \eqn{\mathbb{S}}.
#' @returns \code{data_pattern} A vector where the data are grouped according to \eqn{\mathbb{S}}.
#' @returns \code{muS} The sequence of means.
#' @returns \code{C_S} The sequence of covariance matrices.
#' @returns \code{sigma_squared_S} The sequence of variances.
#' @returns \code{SigmaS} The sequence of correlation matrices.
#' @returns \code{ambient_dimension} The dimension \eqn{d} of the data.
#'
#' @references \insertRef{BB2024}{MCARtest}
#'
#' @importFrom pracma sqrtm
#' @importFrom misty na.indicator
#' @importFrom stats cov2cor var
#'
#' @examples
#' library(copula)
#' library(missMethods)
#' library(misty)
#' n = 1000
#'
#' cp = claytonCopula(param = c(1), dim = 5)
#' P = mvdc(copula = cp, margins = c("exp", "exp", "exp", "exp", "exp"),
#'          paramMargins = list(list(1), list(1), list(1), list(1), list(1)))
#' X = rMvdc(n, P)
#' X = delete_MCAR(X, 0.1, c(1,4,5))
#'
#' get_SigmaS(X)

get_SigmaS = function(X){

  #### create vector with indicators of NA's
  d = dim(X)[2]
  v = 1:d
  pattern_indicator = unique(na.indicator(X))
  n_pattern = dim(pattern_indicator)[1]

  ### create sequence of patterns e.g. c(2,3) if d=3 and just X1 is missing
  patterns = list()
  for (row in 1:n_pattern){
    patterns[[row]] = v[as.logical(pattern_indicator[row,])]
  }

  ### add a column to X, corresponding to the pattern of missingness
  tmp = na.indicator(X)
  extra_col = numeric(length = dim(tmp)[1])
  for (row in 1:dim(tmp)[1]){
    extra_col[row] = paste(v[as.logical(tmp[row,])], collapse = "")
  }

  X = cbind(X, extra_col)

  ### split the data X into multiple subset, each of which corresponds to a certain pattern
  data_pattern = list()
  for (i in 1:n_pattern){
    tmp = matrix(X[X[,dim(X)[2]] == paste(patterns[[i]], collapse = ""),], ncol = dim(X)[2])
    data_pattern[[i]] = apply(matrix(tmp[, patterns[[i]]], ncol = length(patterns[[i]])), c(1,2), as.numeric)
  }

  deletion = c()
  for (i in 1:n_pattern){
    if (dim(data_pattern[[i]])[1] <= dim(data_pattern[[i]])[2]){ ####### remove pattern with sample size too small
      print("Warning: Some patterns were deleted because the sample size was too small!")
      deletion = c(deletion, i)
    }
  }

  if (length(deletion) > 0){
    data_pattern = data_pattern[-deletion]
    patterns = patterns[-deletion]
  }

  n_pattern = n_pattern - length(deletion)

  muS = list()
  C_S = list()
  sigma_squared_S = list()
  SigmaS = list()
  for (i in 1:n_pattern){
    muS[[i]] = colMeans(data_pattern[[i]])
    C_S[[i]] = var(data_pattern[[i]])
    sigma_squared_S[[i]] = diag(C_S[[i]])
    SigmaS[[i]] = cov2cor(C_S[[i]])
  }

  my_list = list("patterns" = patterns, "n_pattern" = n_pattern,
                 "data_pattern" = data_pattern, "muS" = muS,
                 "C_S" = C_S, "sigma_squared_S" = sigma_squared_S,
                 "SigmaS" = SigmaS, "ambient_dimension" = d)
  return(my_list)
}
