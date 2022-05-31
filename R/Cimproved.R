#' Calculate the critical value for our improved test
#'
#' Calculate a critical value for an MCAR test based on knowledge of the facet
#' structure of the Minkowski sum calculated by \code{ConsMinkSumHrep}.
#'
#' @param nS A vector of sample sizes, with each entry corresponding to an observation pattern.
#' @param bS A binary matrix specifying the set of observation patterns. Each row encodes a single pattern.
#' @param M A vector of positive integers giving the alphabet sizes of the discrete variables.
#' @param DR The quantity \eqn{D_R} appearing in \insertCite{BS2022;textual}{MCARtest}.
#' @param Fp The quantity \eqn{F'} appearing in \insertCite{BS2022;textual}{MCARtest}.
#' @param alpha The desired significance level \eqn{\alpha} of the test.
#'
#' @return The critical value \eqn{C_\alpha'} defined in \insertCite{BS2022;textual}{MCARtest}.
#' @export
#'
#' @references \insertRef{BS2022}{MCARtest}
#'
#' @examples
#' bS=matrix(c(1,1,0, 1,0,1, 0,1,1),byrow=TRUE,ncol=3)
#' r=4; s=3
#' M=c(r,s,2)
#' Cimproved(rep(1000,3),bS,M,1,(2^r-2)*(2^s-2),0.05)

Cimproved=function(nS,bS,M,DR,Fp,alpha){
  cS=nrow(bS)

  part1=2*DR^2*log(max(1,2*Fp*cS/alpha))/min(nS)

  part2=0
  for(s1 in 1:(cS-1)){
    for(s2 in (s1+1):cS){
      n=min(nS[s1],nS[s2])

      S1=bS[s1,]; S2=bS[s2,]
      Scap=S1*S2
      ScapDims=Scap*M; ScapDims[ScapDims==0]=1
      cXScap=prod(ScapDims)

      part2=max(part2, (cXScap*log(2)+log(2*cS*(cS-1)/alpha))/n)
    }
  }
  part2=part2*2^(2*cS+7)

  return(cS*sqrt(max(part1,part2)))
}
