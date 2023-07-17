#' Calculate the H-representation of the general (possibly inconsistent) Minkowski sum
#'
#' Computes the minimal halfspace representation of the Minkowski sum of the marginal polytope and the inconsistent ball
#' defined in \insertCite{BS2022;textual}{MCARtest}.
#'
#' @param bS A binary matrix specifying the set of observation patterns. Each row encodes a single pattern.
#' @param M A vector of positive integers giving the alphabet sizes of the discrete variables.
#' @param round A logical value indicating whether or not to round coefficients to 15 significant figures.
#' The function \code{RoundErrors} can be used separately to substitute other values for 15. Defaults to FALSE.
#'
#' @return A halfspace representation object as used by the \code{rcdd} package. See \insertCite{RCDD;textual}{MCARtest} for more detail.
#' @export
#'
#' @references
#' \insertRef{BS2022}{MCARtest}
#'
#' \insertRef{RCDD}{MCARtest}
#'
#' @importFrom rcdd q2d
#' @importFrom rcdd d2q
#' @importFrom rcdd makeH
#' @importFrom rcdd makeV
#' @importFrom rcdd redundant
#' @importFrom rcdd scdd
#' @importFrom rcdd addVpoints
#'
#' @examples
#' bS=matrix(c(1,1, 1,0),byrow=TRUE,ncol=2)
#' InconsMinkSumHrep(bS,c(2,2))

InconsMinkSumHrep=function(bS,M,round=FALSE){
  A=Amatrix(bS,M); cardS=nrow(bS)
  MPolyHrep=MargPolyHrep(bS,M)
  v=ncol(MPolyHrep)-2

  # The set of all (possibly inconsistent) margins
  SumToOneConstr=matrix(rep(0,cardS*v),nrow=cardS); current=1
  for(S in 1:cardS){
    MS=M[bS[S,]==1]
    SumToOneConstr[S,current:(current+prod(MS)-1)]="1"
    current=current+prod(MS)
  }
  PHrep=makeH(a1=-diag(v),b1=rep(0,v),a2=SumToOneConstr,b2=rep(1,cardS))
  PVrep=scdd(PHrep)$output
  PVrep=addVpoints(rep(0,v),PVrep)


  # Combining everything for P_S^{0,*}+P_S^{**}
  Vrep=d2q(makeV(points=q2d(PVrep[,-c(1:2)]), rays=t(A)))
  Vrep=redundant(Vrep)$output

  Hrep=scdd(Vrep)$output
  Hrep=redundant(Hrep)$output
  colnames(Hrep)=colnames(MPolyHrep)

  if(round==FALSE){
    return(Hrep)
  }
  if(round==TRUE){
    return(RoundErrors(Hrep))
  }
}
