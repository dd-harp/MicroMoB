/*
 * sampling random variates
 * Sean L. Wu (@slwu89)
 * 2021-1-11
 */

#include "utils_sample.h"

static inline double get_beta_1_b(const double b){return 1.0 - pow(unif_rand(), 1.0/b);}

SEXP C_draw_multinom(
    SEXP n,
    SEXP probs
) {
  int nn = Rf_asInteger(n);
  double* probs_p = REAL(probs);
  int n_bins = Rf_length(probs);
  SEXP out = Rf_protect(Rf_allocVector(INTSXP, n_bins));
  int* out_p = INTEGER(out);

  GetRNGstate();
  draw_multinom_internal(nn, probs_p, n_bins, out_p);

  PutRNGstate();
  Rf_unprotect(1);
  return out;
}


// this algorithm is from: Startek, Michał. "An asymptotically optimal, online algorithm for weighted random sampling with replacement." arXiv preprint arXiv:1611.00532 (2016).
void draw_multinom_internal(
    int n,
    double* probs,
    int n_bins,
    int* result
) {

  static const double switchover = 1.0;

  for (int i = 0; i < n_bins; ++i) {
    result[i] = 0;
  }

  double pprob = 0.0;
  double cprob = 0.0;
  int pidx = 0;
  while (n > 0) {
    pprob += probs[pidx];
    while (((pprob - cprob) * n / (1.0 - cprob)) < switchover) {
      cprob += get_beta_1_b(n) * (1.0 - cprob);
      while (pprob < cprob) {
        pprob += probs[++pidx];
      }
      if (n_bins == pidx) {
        result[pidx] = 1;
      } else {
        result[pidx] += 1;
      }
      n--;
      if (n == 0) {
        break;
      }
    }
    if (n == 0) {
      break;
    }
    double p = (pprob-cprob)/(1.0-cprob);
    int nrtaken = 0;
    if (p >= 1.0) {
      nrtaken = n;
    } else {
      if (p > 0.0) {
        nrtaken = (int)Rf_rbinom((double)n, p);
      }
    }
    if (nrtaken > 0) {
      if (n_bins == pidx) {
        result[pidx] = nrtaken;
      } else {
        result[pidx] += nrtaken;
      }
    }
    n -= nrtaken;
    pidx++;
    cprob = pprob;
  }

}

