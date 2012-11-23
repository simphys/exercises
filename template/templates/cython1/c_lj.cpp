#include <math.h>
#include <iostream>
using namespace std;

extern "C" {
  int N;
  double L;
  double rcut;
  double shift;

  // this function is required so that one can set the global
  // variables required by the other functions via the Python code
  void c_set_globals(double _L, int _N, double _rcut, double _shift) {
    N = _N;
    L = _L;
    rcut = _rcut;
    shift = _shift;
  }

  // compute the minimum image of particle i and j
  double minimum_image(double *x, int i, int j, double rij[3]) {
    rij[0] = x[j] - x[i];
    rij[1] = x[j+N] - x[i+N];
    rij[2] = x[j+2*N] - x[i+2*N];

    rij[0] -= rint(rij[0]/L)*L;
    rij[1] -= rint(rij[1]/L)*L;
    rij[2] -= rint(rij[2]/L)*L;
  }

  void compute_lj_force(double rij[3], double fij[3]) {
    double r2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
    if (r2 < rcut*rcut) {
      double fac = 4.0 * (12.0 * pow(r2, -6.5) - 6.0 * pow(r2, -3.5)) / sqrt(r2);
      fij[0] = fac*rij[0];
      fij[1] = fac*rij[1];
      fij[2] = fac*rij[2];
    } else {
      fij[0] = 0.0;
      fij[1] = 0.0;
      fij[2] = 0.0;
    }
  }
  
  double compute_lj_potential(double rij[3]) {
    double r2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
    if (r2 < rcut*rcut) {
      return 4.0 * (pow(r2, -6) - pow(r2, -3)) - shift;
    } else {
      return 0.0;
    }
  }

  // this function was prefixed with "c_" to avoid a nameclash with
  // the corresponding Python function
  void c_compute_forces(double *x, double *f) {
    double rij[3];
    double fij[3];

    // set all forces to zero
    for (int i = 0; i < 3*N; i++)
      f[i] = 0.0;

    // add up the forces
    for (int i = 1; i < N; i++) {
      for (int j = 0; j < i; j++) {
        minimum_image(x, i, j, rij);
        compute_lj_force(rij, fij);

        f[i] -= fij[0];
        f[i+N] -= fij[1];
        f[i+2*N] -= fij[2];

        f[j] += fij[0];
        f[j+N] += fij[1];
        f[j+2*N] += fij[2];
      }
    }
  }

  /*
  double c_compute_energy(double *x, double *v) {
  // TODO
  }
  */

}
