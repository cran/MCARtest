#' Generate the matrix A, whose columns are the vertices of the marginal polytope.
#'
#' @param bS A binary matrix specifying the set of observation patterns. Each row encodes a single pattern.
#' @param M A vector of positive integers giving the alphabet sizes of the discrete variables.
#'
#' @return The matrix A.
#' @export
#'
#' @examples
#' bS=matrix(c(1,1,0, 1,0,1, 0,1,1),byrow=TRUE,ncol=3)
#' M=c(2,2,2)
#' Amatrix(bS,M)

Amatrix=function(bS,M){
  X=prod(M)  # Dimension of cX

  bSwDims=t(t(bS)*M); bSwDims[bSwDims==0]=1
  v=sum(apply(bSwDims,1,prod)) # Dimension of cX_bS

  A=matrix(rep(0,X*v),ncol=X)

  cardS=nrow(bS)
  for(i in 1:prod(M)){
    x=as.vector(arrayInd(i,M))
    for(S in 1:cardS){
      xS=x[bS[S,]==1]
      A[row_index(bS,M,S,xS),col_index(M,x)]=1
    }
  }
  return(A)
}
