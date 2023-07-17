#' Round errors in halfspace representations
#'
#' @param X A halfspace representation to be rounded.
#' @param digits An integer giving the number of significant figures to be kept.
#'
#' @return A rounded halfspace representation.
#' @export
#'
#' @importFrom rcdd q2d
#' @importFrom rcdd d2q
#'
#' @examples
#' bS=matrix(c(1,1,1,0, 1,0,0,1, 0,1,0,1, 0,0,1,1),byrow=TRUE,ncol=4)
#' RoundErrors("9007199254740992/6004799503160661") #Occurs in ConsMinkSumHrep(bS,c(2,2,2,2))
#'
#'
#'

RoundErrors=function(X,digits=15){
  return(d2q(signif(q2d(X),digits)))
}
