#include <Rcpp.h>
#include <vector>
#include <set>

using namespace Rcpp;

//' Generates the row indices used internally to generate the sparse matrix A.
//'
//' @param bS A binary matrix specifying the set of observation patterns. Each row encodes a single pattern.
//' @param M A vector of positive integers giving the alphabet sizes of the discrete variables.
//' @return A vector of row indices.
// [[Rcpp::export]]
std::vector<unsigned long long> aMatrixSparseRevLex(std::vector<int>& bS, std::vector<int>& M) {
  int d = M.size();
  unsigned long long cardS = bS.size() / d;

  std::vector<unsigned long long> partialProdRev(d, 1);
  for (int i = 1; i < d; i++) {
    partialProdRev[i] = partialProdRev[i - 1] * M[i - 1];
  }
  unsigned long long cardChi = M[d - 1] * partialProdRev[d - 1];

  std::vector<std::vector<unsigned long long>> partProdRevS(cardS, std::vector<unsigned long long>{});
  std::vector<unsigned long long> szS(cardS, 1);
  for (unsigned long long i = 0; i < cardS; i++) {
    std::vector<int> mS;
    for (int j = 0; j < d; j++) {
      if (bS[i * d + j]) {
        mS.push_back(M[j]);
      }
    }
    partProdRevS[i].resize(mS.size());
    partProdRevS[i][0] = 1;
    for (int j = 1; j < mS.size(); j++) {
      partProdRevS[i][j] = partProdRevS[i][j - 1] * mS[j - 1];
    }
    szS[i] = partProdRevS[i][mS.size() - 1] * mS[mS.size() - 1];
  }

  std::vector<unsigned long long> A(cardChi * cardS);

  for (unsigned long long iChi = 0; iChi < cardChi; iChi++) {
    std::vector<unsigned long long> vChi(d, 0);
    for (unsigned long long i = d, iChiCC = iChi; i > 0; i--) {
      vChi[i - 1] = iChiCC / partialProdRev[i - 1];
      iChiCC -= vChi[i - 1] * partialProdRev[i - 1];
    }

    unsigned long long partSz = 0;
    for (unsigned long long i = 0; i < cardS; i++) {
      unsigned long long indS = 0;
      for (int j = d - 1, k = (int) partProdRevS[i].size() - 1; j >= 0; j--) {
        if (bS[i * d + j]) {
          indS += vChi[j] * partProdRevS[i][k];
          k--;
        }
      }
      A[iChi * cardS + i] = partSz + indS + 1;
      partSz += szS[i];
    }
  }

  return A;
}

//' Calculates the total cardinality of the sample spaces.
//'
//' @param bS A binary matrix specifying the set of observation patterns. Each row encodes a single pattern.
//' @param M A vector of positive integers giving the alphabet sizes of the discrete variables.
//' @return The total cardinality.
// [[Rcpp::export]]
unsigned long long infoS(std::vector<int>& bS, std::vector<int>& M) {
  unsigned long long ans = 0;
  int d = M.size();
  unsigned long long cardS = bS.size() / d;

  for (unsigned long long i = 0; i < cardS; i++) {
    unsigned long long szs = 1;
    for (int j = 0; j < d; j++) {
      if (bS[i * d + j]) {
        szs *= M[j];
      }
    }
    ans += szs;
  }

  return ans;
}

//' Calculates the individual cardinalities of the sample spaces.
//'
//' @param bS A binary matrix specifying the set of observation patterns. Each row encodes a single pattern.
//' @param M A vector of positive integers giving the alphabet sizes of the discrete variables.
//' @return A vector of individual cardinalities.
// [[Rcpp::export]]
std::vector<unsigned long long> infoS2(std::vector<int>& bS, std::vector<int>& M) {
  int d = M.size();
  unsigned long long cardS = bS.size() / d;
  std::vector<unsigned long long> szs(cardS, 1);

  for (unsigned long long i = 0; i < cardS; i++) {
    for (int j = 0; j < d; j++) {
      if (bS[i * d + j]) {
        szs[i] *= M[j];
      }
    }
  }

  return szs;
}

//' Generates the column indices used internally to generate the sparse matrix A.
//'
//' @param cardS The number of missingness patterns.
//' @param cardChi The cardinality of the full joint space.
//' @return A vector of column indices.
// [[Rcpp::export]]
std::vector<unsigned long long> colVector(unsigned long long cardS, unsigned long long cardChi) {
  std::vector<unsigned long long> v(cardS * cardChi, 0);
  for (unsigned long long j = 0; j < cardChi; j++) {
    for (unsigned long long i = 0; i < cardS; i++) {
      v[j * cardS + i] = j + 1;
    }
  }
  return v;
}

//' Internal function multiplying a mass function by the sparse matrix A.
//' 
//' @param p A subprobability mass function on the full joint space.
//' @param bS A binary matrix specifying the set of observation patterns. Each row encodes a single pattern.
//' @param M A vector of positive integers giving the alphabet sizes of the discrete variables.
//' @return A collection of marginal mass functions.
// [[Rcpp::export]]
std::vector<double> margProj(std::vector<double>& p, std::vector<int>& bS, std::vector<int>& M) {
  unsigned long long cardS = bS.size() / M.size();
  unsigned long long cardChi = p.size();
  auto infobS = infoS(bS, M);
  auto infobS2 = infoS2(bS, M);

  auto A = aMatrixSparseRevLex(bS, M);
  std::vector<std::set<unsigned long long>> ARow(infobS, std::set<unsigned long long>{});
  for (unsigned long long i = 0; i < cardS * cardChi; i++) {
    ARow[A[i] - 1].insert(i / cardS);
  }

  std::vector<double> projpS(infobS, 0);
  for (unsigned long long i = 0; i < infobS; i++) {
    for (auto j : ARow[i]) {
      projpS[i] += p[j];
    }
  }

  return projpS;
}


