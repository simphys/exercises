#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

int N;
double L;
double rcut;
double shift;

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

void compute_forces(double *x, double *f) {
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
double compute_energy(double *x, double *v) {
  // TODO:
  // .
  // .
  // .
}
*/

void step_vv(double *x, double *v, double *f, double dt) {
  for (int i = 0; i < N; i++) {
    x[i] += v[i]*dt + 0.5*f[i] * dt*dt;
    x[i+N] += v[i+N]*dt + 0.5*f[i+N] * dt*dt;
    x[i+2*N] += v[i+2*N]*dt + 0.5*f[i+2*N] * dt*dt;
  }

  for (int i = 0; i < N; i++) {
    v[i] += 0.5*f[i] * dt;
    v[i+N] += 0.5*f[i+N] * dt;
    v[i+2*N] += 0.5*f[i+2*N] * dt;
  }

  compute_forces(x, f);

  for (int i = 0; i < N; i++) {
    v[i] += 0.5*f[i] * dt;
    v[i+N] += 0.5*f[i+N] * dt;
    v[i+2*N] += 0.5*f[i+2*N] * dt;
  }
}

int main() {
  // Constants
  const double density = 0.7;
  const int n = 3;
  const double dt = 0.01;
  const double tmax = 1.0;
  rcut = 2.5;
  shift = -0.016316891136;

  N = n*n*n;
  const double volume = N/density;
  L = pow(volume , 1./3.);
  const double l = L/n;

  double x[3*N];
  double v[3*N];
  double f[3*N];

  double t = 0.0;

  /* Set up cubic lattice */
  for (int i = 0; i < n; i++) {
    double x0 = i*l;
    for (int j = 0; j < n; j++) {
      double x1 = j*l;
      for (int k = 0; k < n; k++) {
        double x2 = k*l;
        int ix = i*n*n + j*n + k;
        x[ix] = x0;
        x[ix+N] = x1;
        x[ix+2*N] = x2;

        v[ix] = drand48()*2.0 - 1.0;
        v[ix+N] = drand48()*2.0 - 1.0;
        v[ix+2*N] = drand48()*2.0 - 1.0;
      }
    }
  }
  
  compute_forces(x, f);

  ofstream vtffile("ljfluid_c.vtf");
  vtffile << "atom 0:" << (N-1) << " radius 0.5" << endl;
  vtffile << "unitcell " << L << " " << L << " " << L << endl;

  // ofstream datfile("ljfluid_c.dat");

  while (t < tmax) {
    step_vv(x, v, f, dt);
    t += dt;

    // TODO: implement compute_energy!
    //    double E = compute_energy(x, v);
    // datfile << t << '\t' << E << endl;
    cout << t << endl;

    vtffile << "timestep" << endl;
    for (int i = 0; i < N; i++) {
      vtffile << x[i] << '\t' 
              << x[i+N] << '\t' 
              << x[i+2*N] << endl;
    }
  }

  // datfile.close();
  vtffile.close();

}
