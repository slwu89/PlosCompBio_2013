#include <R.h>
#include <Rinternals.h>

#include <stdio.h>

#include <gsl/gsl_rng.h>

SEXP hello_C(SEXP seedR){
  unsigned long int seed = asInteger(seedR);
  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,seed);
  double aa = gsl_rng_uniform_pos(r);
  printf("hello world! here is a random number from the gsl ... %g\n",aa);
  return R_NilValue;
};
