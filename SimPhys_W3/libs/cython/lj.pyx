import cython

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef struct MyEnergies:
    double E
    double E_pot
    double E_kin
cdef extern void c_set_globals(double _L, int _N, double _rcut, double _shift)
cdef extern void c_compute_forces(double *x, double *f)
cdef extern MyEnergies c_compute_energy (double *x, double *v)
cdef extern double c_rebuild_neighbor_lists(double *x, double vlsize)

def set_globals(double L, int N, double rcut, double shift):
    c_set_globals(L, N, rcut, shift)

def compute_forces(np.ndarray[double, ndim=2, mode='c'] x not None):
    cdef np.ndarray[double, ndim=2, mode='c'] f = np.zeros_like(x)
    N = x.shape[1]
    c_compute_forces (&x[0,0], &f[0,0])
    return f

def compute_energy(np.ndarray[double, ndim=2, mode='c'] x not None, 
                   np.ndarray[double, ndim=2, mode='c'] v not None):
    N = x.shape[1]
    r = c_compute_energy(&x[0,0], &v[0,0])
    return r.E, r.E_pot, r.E_kin

def rebuild_neighbor_lists(np.ndarray[double, ndim=2, mode='c'] x not None,
                           double vlsize):
    N = x.shape[1]
    c_rebuild_neighbor_lists(&x[0,0], vlsize)
