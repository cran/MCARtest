#' Generate the matrix A, whose columns are the vertices of the marginal polytope, as a sparse matrix.
#'
#' @param bS A binary matrix specifying the set of observation patterns. Each row encodes a single pattern.
#' @param M A vector of positive integers giving the alphabet sizes of the discrete variables.
#'
#' @return The matrix A.
#' @export
#'
#' @importFrom Rcpp evalCpp
#' @importFrom Matrix sparseMatrix
#'
#' @examples
#' bS=matrix(c(1,1,0, 1,0,1, 0,1,1),byrow=TRUE,ncol=3)
#' M=c(2,2,2)
#' AmatrixSparse(bS,M)
#'
#'
AmatrixSparse <- function(bS, M) {
  cardS <- nrow(bS)
  cardChi <- prod(M)

  bS <- c(t(bS))
  ii <- aMatrixSparseRevLex(bS, M)
  jj <- colVector(cardS, cardChi)
  totCardS <- infoS(bS, M)
  a <- Matrix::sparseMatrix(i = ii, j = jj, x = 1, dims = c(totCardS, cardChi))

  return(a)
}
