#' Carry out Little's test of MCAR
#'
#' @param X The dataset with incomplete data, where all the pairs of variables are observed together.
#' @param alpha The nominal level of the test.
#' @param type Determines the test statistic to use, based on the discussion in Section 5 in \insertCite{BB2024;textual}{MCARtest}.
#' The default option is "mean&cov", and uses the test statistic \eqn{d^2_{\mathrm{aug}}}. When set equal to "cov", implements a test
#' of MCAR based on \eqn{d^2_{\mathrm{cov}}}, while, when set equal to "mean", implements the classical Little's test as defined in
#' \insertCite{Little1988;textual}{MCARtest}.
#'
#' @references \insertRef{BB2024}{MCARtest}
#' @references \insertRef{Little1988}{MCARtest}
#'
#' @return A Boolean, where TRUE stands for reject MCAR. This is computed by comparing the p-value of Little's test,
#' found by comparing the log likelihood ratio statistic to the chi-squared distribution with the appropriate number
#' of degrees of freedom, with the nominal level \code{alpha}. Described in \insertCite{Little1988;textual}{MCARtest}.
#' @export
#'
#' @importFrom norm prelim.norm em.norm getparam.norm
#' @importFrom stats cor qchisq
#'
#' @examples
#' library(MASS)
#' alpha = 0.05
#' n = 200
#'
#' SigmaS=list() #Random 2x2 correlation matrices (necessarily consistent)
#' for(j in 1:3){
#' x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1)
#' SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
#' }
#'
#' X1 = mvrnorm(n, c(0,0), SigmaS[[1]])
#' X2 = mvrnorm(n, c(0,0), SigmaS[[2]])
#' X3 = mvrnorm(n, c(0,0), SigmaS[[3]])
#' columns = c("X1","X2","X3")
#' X = data.frame(matrix(nrow = 3*n, ncol = 3))
#' X[1:n, c("X1", "X2")] = X1
#' X[(n+1):(2*n), c("X2", "X3")] = X2
#' X[(2*n+1):(3*n), c("X1", "X3")] = X3
#' X = as.matrix(X)
#'
#' little_test(X, alpha)

little_test = function(X, alpha, type="mean&cov"){
  s = prelim.norm(as.matrix(X))
  thetahat = em.norm(s)
  mu_true = getparam.norm(s,thetahat,corr=TRUE)$mu; Sigma_true = getparam.norm(s,thetahat,corr=TRUE)$r

  result = get_SigmaS(X)
  SigmaS = result$SigmaS; patterns = result$pattern; n_pattern = result$n_pattern; data_pattern = result$data_pattern
  d = result$ambient_dimension

  if(type == "mean"){
    d_squared = 0
    df = -d

    for (i in 1:n_pattern){
      n_S = dim(data_pattern[[i]])[1]
      card_S = dim(data_pattern[[i]])[2]

      x_S = colMeans(data_pattern[[i]]) - mu_true[patterns[[i]]]
      L_S = Sigma_true[patterns[[i]],patterns[[i]]]
      Sigma_S = cor(data_pattern[[i]])

      d_squared = d_squared + n_S*t(x_S)%*%solve(n_S*L_S/(n_S-1))%*%t(t(x_S))
      df = df + card_S
      print(df)
    }
    little_d = (d_squared > qchisq(1-alpha, df))
  }

  else if(type == "cov"){

    for (i in length(SigmaS)){
      if (min(eigen(SigmaS[[i]])$values) < 10^-7){
        print("SigmaS is singular!")
        return(NA)
      }
    }

    d_cov = 0
    df = -d*(d+1)/2

    for (i in 1:n_pattern){
      n_S = dim(data_pattern[[i]])[1]
      card_S = dim(data_pattern[[i]])[2]

      x_S = colMeans(data_pattern[[i]]) - mu_true[patterns[[i]]]
      L_S = Sigma_true[patterns[[i]],patterns[[i]]]
      Sigma_S = cor(data_pattern[[i]])

      d_cov = d_cov + n_S*(sum(diag(Sigma_S%*%solve(L_S))) - card_S - log(abs(det(as.matrix(Sigma_S)))) +
                             log(det(as.matrix(L_S))))
      df = df + card_S*(card_S+1)/2
      print(df)
    }

    little_d = (d_cov > qchisq(1-alpha, df))
  }

  else{

    for (i in length(SigmaS)){
      if (min(eigen(SigmaS[[i]])$values) < 10^-7){
        print("SigmaS is singular!")
        return(NA)
      }
    }

    d_aug = 0
    df = -d*(d+3)/2

    for (i in 1:n_pattern){
      n_S = dim(data_pattern[[i]])[1]
      card_S = dim(data_pattern[[i]])[2]

      x_S = colMeans(data_pattern[[i]]) - mu_true[patterns[[i]]]
      L_S = Sigma_true[patterns[[i]],patterns[[i]]]
      Sigma_S = cor(data_pattern[[i]])

      d_aug = d_aug + n_S*t(x_S)%*%solve(n_S*L_S/(n_S-1))%*%t(t(x_S)) +
        n_S*(sum(diag(Sigma_S%*%solve(L_S))) - card_S - log(abs(det(as.matrix(Sigma_S)))) +
               log(det(as.matrix(L_S))))
      df = df + card_S*(card_S+3)/2
      print(df)
    }

    little_d = (d_aug > qchisq(1-alpha, df))
  }

  return(little_d)

}
