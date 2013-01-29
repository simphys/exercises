import cython
import numpy
cimport numpy

# declare the interface to the C code
cdef extern void c_init()
cdef extern void c_random_list(int N, double *rs)

# Cython interface to C function
def random_list(N):
    # define an (empty) numpy double array to hold the data
    cdef numpy.ndarray[double, ndim=1, mode='c'] rs = numpy.empty((N))
    # pass the array to the C function
    c_random_list(N, &rs[0])
    return rs

# call c_init() when the module is loaded to init the gsl RNG
c_init()
