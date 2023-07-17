#' Calculate the H-representation of the marginal polytope
#'
#' Computes the minimal halfspace representation of the marginal polytope
#' defined, for example, in \insertCite{BS2022;textual}{MCARtest}.
#'
#' @param bS A binary matrix specifying the set of observation patterns. Each row encodes a single pattern.
#' @param M A vector of positive integers giving the alphabet sizes of the discrete variables.
#' @param round A logical value indicating whether or not to round coefficients to 15 significant figures.
#' The function \code{RoundErrors} can be used separately to substitute other values for 15. Defaults to FALSE.
#'
#' @return A halfspace representation object as used by the \code{rcdd} package. See \insertCite{RCDD;textual}{MCARtest} for more detail.
#' @export
#'
#' @importFrom rcdd d2q
#' @importFrom rcdd makeV
#' @importFrom rcdd redundant
#' @importFrom rcdd scdd
#'
#' @references
#' \insertRef{BS2022}{MCARtest}
#'
#' \insertRef{RCDD}{MCARtest}
#'
#' @examples
#' bS=matrix(c(1,1,0, 1,0,1, 0,1,1),byrow=TRUE,ncol=3)
#' MargPolyHrep(bS,c(2,2,2))


MargPolyHrep=function(bS,M,round=FALSE){
  A=Amatrix(bS,M)
  v=nrow(A); cardS=nrow(bS)

  Vrep=d2q(makeV(points=t(A)))
  Hrep=scdd(Vrep)$output
  Hrep=redundant(Hrep)$output

  Names=rep(0,v); Names[1:2]=c("l","b"); current=2
  for(S in 1:cardS){
    MS=M[bS[S,]==1]
    for(i in 1:prod(MS)){
      current=current+1
      xS=as.vector(arrayInd(i,MS))
      Name=bS[S,]; Name[Name==0]="d"; Name[Name==1]=xS
      Name=paste(Name,collapse="")
      Names[current]=paste("p",Name,sep="")
    }
  }
  colnames(Hrep)=Names
  if(round==FALSE){
    return(Hrep)
  }
  if(round==TRUE){
    return(RoundErrors(Hrep))
  }

}
