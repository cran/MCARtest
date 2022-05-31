#' Perform one step of the EM algorithm for finding the MLE under MCAR in a contingency table.
#'
#' @param pt An input probability mass function on the joint space, to be updated.
#' @param p0h An empirical mass function calculated using complete observations.
#' @param n0 An integer giving the number of complete observations used to calculate \code{p0h}.
#' @param pSh A sequence of empirical mass functions calculated using incomplete observations.
#' @param nS A sequence of integers giving the numbers of incomplete observations used to calculate \code{pSh}.
#' @param bS A binary matrix specifying the set of observation patterns. Each row encodes a single pattern.
#' @param M A vector of positive integers giving the alphabet sizes of the discrete variables.
#'
#' @return The updated probability mass function on the joint space.
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
#' EMiteration(p0,p0h,n0,pSh,nS,bS,M)
#'

EMiteration=function(pt,p0h,n0,pSh,nS,bS,M){
  X=p0h*n0
  for(S in 1:nrow(bS)){
    MS=M[bS[S,]==1]
    start=row_index(bS,M,S,rep(1,sum(bS[S,]))); end=row_index(bS,M,S,MS)
    ps=apply(pt,which(bS[S,]==1),sum)
    Xs=array(nS[S]*pSh[start:end],dim=MS)
    xs=outer(Xs/ps, array(1,dim=M[bS[S,]==0]))
    permvec=rep(1,ncol(bS))
    permvec[bS[S,]==1]=(1:length(MS))
    permvec[bS[S,]==0]=((length(MS)+1):length(M))
    xs=aperm(xs,perm=permvec)
    X=X+xs*pt
  }
  return(X/(n0+sum(nS)))
}
