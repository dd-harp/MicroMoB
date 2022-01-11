/*
 * sampling random variates
 * Sean L. Wu (@slwu89)
 * 2021-1-11
 */

#ifndef UTILS_SAMPLE_H
#define UTILS_SAMPLE_H

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

SEXP draw_multinom_(
    SEXP n,
    SEXP probs
);

// n: number of balls
// probs: pointer to weights of bins
// n_bins: number of bins
// result: pointer to array to write result
void draw_multinom_internal(
  int n,
  double* probs,
  int n_bins,
  int* result
);



#endif
