#' Calculate the H-representation of the consistent Minkowski sum
#'
#' Computes the minimal halfspace representation of the Minkowski sum of the marginal polytope
#' and the consistent ball defined in \insertCite{BS2022;textual}{MCARtest}.
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
#' bS=matrix(c(1,1,0, 1,0,1, 0,1,1),byrow=TRUE,ncol=3)
#' ConsMinkSumHrep(bS,c(2,2,2))
#'


ConsMinkSumHrep=function(bS,M,round=FALSE){
  A=Amatrix(bS,M)
  MPolyHrep=MargPolyHrep(bS,M)
  v=ncol(MPolyHrep)-2

  ConsConstr=MPolyHrep[MPolyHrep[,1]=="1",-c(1,2)]
  ConsRHS=d2q(-q2d(MPolyHrep[MPolyHrep[,1]=="1",2]))

  # The set of all consistent submargins (P_S^{cons,**})
  PConsHrep=makeH(a1=-diag(v),b1=rep(0,v),a2=ConsConstr,b2=ConsRHS)
  PConsVrep=scdd(PConsHrep)$output
  PConsVrep=addVpoints(rep(0,v),PConsVrep)


  # Combining everything for P_S^{0,*}+P_S^{cons,**}
  Vrep=d2q(makeV(points=q2d(PConsVrep[,-c(1:2)]), rays=t(A)))  # The V-representation of the marginal cone plus the set of valid margins
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
