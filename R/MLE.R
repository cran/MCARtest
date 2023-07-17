#' Compute the MLE under MCAR in a contingency table using the EM algorithm, given complete and incomplete observations.
#'
#' @param p0h An empirical mass function calculated using complete observations.
#' @param n0 An integer giving the number of complete observations used to calculate \code{p0h}.
#' @param pSh A sequence of empirical mass functions calculated using incomplete observations.
#' @param nS A sequence of integers giving the numbers of incomplete observations used to calculate \code{pSh}.
#' @param bS A binary matrix specifying the set of observation patterns. Each row encodes a single pattern.
#' @param M A vector of positive integers giving the alphabet sizes of the discrete variables.
#' @param Niter An integer giving the number of iterations to be used in the EM algorithm.
#' @param loglik A logical value indicating whether or not the log likelihoods at each step of the EM algorithm should be an output. Defaults to FALSE.
#'
#' @return The output of the EM algorithm, approximating the MLE for the probability mass function on the joint space.
#' @export
#'
#' @importFrom stats runif
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
#' MLE(p0h,n0,pSh,nS,bS,M,50)
#'
#' trace=MLE(p0h,n0,pSh,nS,bS,M,50,loglik=TRUE)[[2]]
#' plot(1:50,trace,type="l")

MLE=function(p0h,n0,pSh,nS,bS,M,Niter,loglik=FALSE){
  pt=array(runif(prod(M)),dim=M); pt=pt/sum(pt)
  if(loglik==TRUE) ll=rep(0,Niter)
  for(niter in 1:Niter){
    pt=EMiteration(pt,p0h,n0,pSh,nS,bS,M)
    if(loglik==TRUE) ll[niter]=loglik0(pt,p0h,n0,pSh,nS,bS,M)
  }
  if(loglik==FALSE){
    return(pt)
  }else{
    return(list(pt,ll))
  }
}
