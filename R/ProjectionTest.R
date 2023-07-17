#' Carry out a test of MCAR in a contingency table, given incomplete observations.
#'
#' @param pSh A sequence of empirical mass functions calculated using incomplete observations.
#' @param nS A sequence of integers giving the numbers of incomplete observations used to calculate \code{pSh}.
#' @param bS A binary matrix specifying the set of observation patterns. Each row encodes a single pattern.
#' @param M A vector of positive integers giving the alphabet sizes of the discrete variables.
#' @param B An integer giving the number of bootstrap samples to be used to calibrate the test.
#'
#' @return The p-value the Monte Carlo test described in \insertCite{BS2022;textual}{MCARtest}.
#' @return The value of the test statistic \eqn{R()}.
#' @export
#'
#' @references
#' \insertRef{BS2022}{MCARtest}
#'
#' @importFrom stats rmultinom
#' @importFrom Rcpp evalCpp
#'
#' @examples
#' bS=matrix(c(1,1,0, 1,0,1, 0,1,1),byrow=TRUE,ncol=3) # Our canonical 3d example
#' M=c(2,2,2)
#' nS=c(200,200,200)
#'
#' pS=c(0.125,0.375,0.375,0.125,0.250,0.250,0.250,0.250,0.100,0.400,0.400,0.100)
#' P12=pS[1:4]; P13=pS[5:8]; P23=pS[9:12]
#' X12=t(rmultinom(1,size=nS[1],prob=P12)/nS[1])
#' X13=t(rmultinom(1,size=nS[2],prob=P13)/nS[2])
#' X23=t(rmultinom(1,size=nS[3],prob=P23)/nS[3])
#' pSh=cbind(X12,X13,X23)
#'
#' ProjectionTest(pSh,nS,bS,M,99)

ProjectionTest <- function(pSh, nS, bS, M, B) {
  TestProg = RindexDual(pSh, bS, M)
  TestStat = TestProg[[1]]
  Projp = TestProg[[2]]

  ProjpS <- margProj(Projp, c(t(bS)), M) / (1 - TestStat)
  szs <- infoS2(c(t(bS)), M)

  stat <- 1;
  for (i in 1:B) {
    boot = numeric(length(pSh))
    cur = 1;
    for (j in 1:nrow(bS)) {
      probS = ProjpS[cur:(cur+szs[j]-1)]
      boot[cur:(cur+szs[j]-1)] = t(rmultinom(1, size = nS[j], prob = probS)/nS[j])
      cur = cur + szs[j]
    }
    bootStat = RindexDual(boot, bS, M)
    if (bootStat[[1]] >= TestStat) {
      stat <- stat + 1
    }
  }

  return(list(stat / (B + 1), TestStat))
}
