#' A function indexing the columns of A
#'
#' A map from the joint space to an index set.
#'
#' @param M A vector of positive integers giving the alphabet sizes of the discrete variables.
#' @param x An element of the joint space.
#'
#' @return A positive integer no greater than the cardinality of the joint space uniquely identifying \code{x}.
#' @export
#'
#' @examples
#' M=c(2,2,2)
#' col_index(M,c(1,1,1))
#' col_index(M,c(1,1,2))
#'
#' M=c(4,3,2)
#' col_index(M,c(1,1,1))
#' col_index(M,c(2,1,1))
#' col_index(M,c(1,2,1))
#' col_index(M,c(1,1,2))

col_index=function(M,x){
  d=length(M)
  Prods=round(exp(lower.tri(diag(d))%*%log(M)))
  return(1+sum((x-1)*Prods))
}
