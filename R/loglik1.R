#' Compute the log likelihood of a probability mass function, without assuming MCAR, given complete and incomplete data
#'
#' @param p0 A probability mass function on the joint space.
#' @param pS A sequence of probability mass functions on the marginal spaces.
#' @param p0h An empirical mass function calculated using complete observations.
#' @param n0 An integer giving the number of complete observations used to calculate \code{p0h}.
#' @param pSh A sequence of empirical mass functions calculated using incomplete observations.
#' @param nS A sequence of integers giving the numbers of incomplete observations used to calculate \code{pSh}.
#' @param bS A binary matrix specifying the set of observation patterns. Each row encodes a single pattern.
#' @param M A vector of positive integers giving the alphabet sizes of the discrete variables.
#'
#' @return The value of the log likelihood.
#' @export
#'
#' @examples
#' bS=matrix(c(1,1,0, 1,0,1, 0,1,1),byrow=TRUE,ncol=3) # Our canonical 3d example
#' M=c(2,2,2)
#' n0=200
#' nS=c(200,200,200)
#'
#' pS=c(0.125,0.375,0.375,0.125,0.250,0.250,0.250,0.250,0.100,0.400,0.400,0.100)
#' P12=pS[1:4]; P13=pS[5:8]; P23=pS[9:12]
#' X12=t(rmultinom(1,size=nS[1],prob=P12)/nS[1])
#' X13=t(rmultinom(1,size=nS[2],prob=P13)/nS[2])
#' X23=t(rmultinom(1,size=nS[3],prob=P23)/nS[3])
#' pSh=cbind(X12,X13,X23)
#'
#' p0=array(0.125,dim=c(2,2,2))
#' p0h=array(rmultinom(1,n0,p0),dim=M)/n0
#'
#' loglik1(p0,pS,p0h,n0,pSh,nS,bS,M)
#'

loglik1=function(p0,pS,p0h,n0,pSh,nS,bS,M){
  ind0=(p0h>0); lik=n0*sum(p0h[ind0]*log(p0[ind0]))
  for(S in 1:nrow(bS)){
    MS=M[bS[S,]==1]
    start=row_index(bS,M,S,rep(1,sum(bS[S,]))); end=row_index(bS,M,S,MS)
    ps=pS[start:end]
    Xs=nS[S]*pSh[start:end]
    inds=(Xs>0)
    lik=lik+sum(Xs[inds]*log(ps[inds]))
  }
  return(lik)
}
