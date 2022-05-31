#' A function indexing the rows of A
#'
#' A map from the observation space to an index set.
#'
#' @param bS A binary matrix specifying the set of observation patterns. Each row encodes a single pattern.
#' @param M A vector of positive integers giving the alphabet sizes of the discrete variables.
#' @param S An integer indicating which observation pattern is of interest.
#' @param xS An element of the observation space of the specified observation pattern.
#'
#' @return A positive integer no larger than the cardinality of the joint space uniquely identifying \code{x}.
#' @export
#'
#' @examples
#' bS=matrix(c(1,1,0, 1,0,1, 0,1,1),byrow=TRUE,ncol=3)
#' M=c(2,2,2)
#' row_index(bS,M,1,c(1,1))
#' row_index(bS,M,2,c(1,1))
#' row_index(bS,M,3,c(1,1))
#'
#' bS=matrix(c(1,1,0, 1,0,1, 0,1,1),byrow=TRUE,ncol=3)
#' M=c(4,3,2)
#' row_index(bS,M,1,c(1,1))
#' row_index(bS,M,1,c(2,1))
#' row_index(bS,M,1,c(3,1))
#' row_index(bS,M,1,c(4,1))
#' row_index(bS,M,1,c(1,2))
#' row_index(bS,M,1,c(2,2))
#'

row_index=function(bS,M,S,xS){
  current=0
  bSwDims=t(t(bS)*M); bSwDims[bSwDims==0]=1
  if(S>1) current=current+sum(apply(bSwDims,1,prod)[1:(S-1)])
  current=current+col_index(M[bS[S,]==1],xS)
  return(current)
}
