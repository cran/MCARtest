#' Calculate the critical value for our simple test
#'
#' Calculate a simple critical value for an MCAR test using only knowledge of
#' the set of observation patterns and the joint observation space.
#'
#' @param nS A vector of sample sizes, with each entry corresponding to an observation pattern.
#' @param bS A binary matrix specifying the set of observation patterns. Each row encodes a single pattern.
#' @param M A vector of positive integers giving the alphabet sizes of the discrete variables.
#' @param alpha The desired significance level \eqn{\alpha} of the test.
#'
#' @return The universal critical value defined in \insertCite{BS2022;textual}{MCARtest}.
#' @export
#'
#' @references \insertRef{BS2022}{MCARtest}
#'
#' @examples
#' bS=matrix(c(1,1,0, 1,0,1, 0,1,1),byrow=TRUE,ncol=3)
#' r=4; s=3
#' M=c(r,s,2)
#' Csimple(rep(1000,3),bS,M,0.05)


Csimple=function(nS,bS,M,alpha){
  bSwDims=t(t(bS)*M); bSwDims[bSwDims==0]=1
  cXS=apply(bSwDims,1,prod) # Dimensions of cX_bS
  return(0.5*sum(sqrt((cXS-1)/nS)) + sqrt(0.5*log(1/alpha)*sum(1/nS)))
}
