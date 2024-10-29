#' Carry out a test of MCAR using first and second moments.
#'
#' This is the implementation of Algorithm 1 in \insertCite{BB2024;textual}{MCARtest}.
#'
#' @param X The dataset with incomplete data.
#' @param alpha The nominal level of the test.
#' @param B The bootstrap sample \eqn{B} for the bootstrap test.
#'
#' @return A Boolean, where TRUE stands for reject MCAR. This is found as outlined in
#' Section 5.2 in \insertCite{BB2024;textual}{MCARtest}.
#' @export
#'
#' @references \insertRef{BB2024}{MCARtest}
#'
#' @importFrom pracma sqrtm
#'
#' @examples
#' library(MASS)
#' alpha = 0.05
#' B = 20
#' m = 500
#' 
#' SigmaS=list() #Random 2x2 correlation matrices (necessarily consistent)
#' for(j in 1:3){
#' x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1)
#' SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
#' }
#' 
#' X1 = mvrnorm(m, c(0,0), SigmaS[[1]])
#' X2 = mvrnorm(m, c(0,0), SigmaS[[2]])
#' X3 = mvrnorm(m, c(0,0), SigmaS[[3]])
#' columns = c("X1","X2","X3")
#' X = data.frame(matrix(nrow = 3*m, ncol = 3))
#' X[1:m, c("X1", "X2")] = X1
#' X[(m+1):(2*m), c("X2", "X3")] = X2
#' X[(2*m+1):(3*m), c("X1", "X3")] = X3
#' X = as.matrix(X)
#' 
#' MCAR_meancovTest(X, alpha, B)

MCAR_meancovTest = function(X, alpha, B){

  n = dim(X)[1]
  result = get_SigmaS(X); d = result$ambient_dimension; av_sigma = compute_av("var", X)

  for (j in 1:d){
    X[,j] = X[,j]/sqrt(av_sigma[j])
  }

  result = get_SigmaS(X)
  av_mu = compute_av("mean", X)
  muS = result$muS; C_S = result$C_S; sigma_squared_S = result$sigma_squared_S; SigmaS = result$SigmaS
  patterns = result$pattern; n_pattern = result$n_pattern; data_pattern = result$data_pattern

  for (i in length(SigmaS)){
    if (min(eigen(SigmaS[[i]])$values) < 10^-7){
      print("SigmaS is singular!")
      return(NA)
    }
  }

  tmp = computeR(patterns, SigmaS)
  T_hat_0 = tmp$R + M(muS, patterns) + V(sigma_squared_S, patterns)
  Q_hat = tmp$Sigma/(1-tmp$R)
  QS_hat = list()
  for (i in 1:n_pattern){
    QS_hat[[i]] = Q_hat[patterns[[i]], patterns[[i]]]
  }

    #### rotate X, to make it look like it's from H0
    rot_data_pattern = list()
    for (i in 1:n_pattern){
      rot_data_pattern[[i]] = t((sqrtm(as.matrix(QS_hat[[i]]))$B)%*%(solve(sqrtm(as.matrix(C_S[[i]]))$B))%*%
                                  t(data_pattern[[i]] - muS[[i]] +  av_mu[patterns[[i]]]))
    }

    sum_indicator = 0
    for (b in 1:B){

      r_ind = 0
      X = data.frame(matrix(nrow = n, ncol = d))
      for (i in 1:n_pattern){
        n_S = dim(rot_data_pattern[[i]])[1]
        tmp_data = as.matrix(rot_data_pattern[[i]][sample(1:n_S, n_S, replace = T),])
        if(dim(unique(tmp_data))[1] < 2){
          tmp_data = as.matrix(rot_data_pattern[[i]][sample(1:n_S, n_S, replace = T),])
        }
        X[(1+r_ind):(n_S+r_ind), patterns[[i]]] = tmp_data
        r_ind = r_ind + n_S
      }
      X = as.matrix(X[1:r_ind,])

      result = get_SigmaS(X); av_sigma = compute_av("var", X)

      for (j in 1:d){
        X[,j] = X[,j]/sqrt(av_sigma[j])
      }

      result = get_SigmaS(X)
      patterns = result$pattern; n_pattern = result$n_pattern; data_pattern = result$data_pattern
      muS_b = result$muS; C_S_b = result$C_S; sigma_squared_S_b = result$sigma_squared_S; SigmaS_b = result$SigmaS

      T_hat_b = computeR(patterns, SigmaS_b)$R + M(muS_b, patterns) + V(sigma_squared_S_b, patterns)

      if (T_hat_b >= T_hat_0){
        sum_indicator = sum_indicator + 1
      }
    }

  p_hat = (1+sum_indicator)/(B+1)
  decision = p_hat < alpha
  return(decision)
}
