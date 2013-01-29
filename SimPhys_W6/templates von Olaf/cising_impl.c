#include <gsl/gsl_rng.h>
#include <stdlib.h>

gsl_rng *r = NULL;

/* Initialize the random number generator. This function is called from cising.pyx. */
void c_init() {
  r = gsl_rng_alloc(gsl_rng_taus);
}

/* Generate a random integer between 0 and N-1. */
int gsl_randint(int N) {
  return gsl_rng_uniform_int(r, N);
}

/* Generate a random float between 0.0 and 1.0. */
double gsl_rand() {
  return gsl_rng_uniform(r);
}

/* This is an example for a function that is exported to Cython and
   uses the GSL RNG. */
void c_random_list(int N, double* rs) {
  int i;
  for (i = 0; i < N; i++)
    rs[i] = gsl_rand();
}
