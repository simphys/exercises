#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
import os, sys, pickle
import numpy as np
import scipy as sp
import scipy.signal
from libs.simlib import Plotter

"""==== PARAMETERS ===="""
# path to simulation File
datafilename = "./data/series.dat"
p = Plotter(show = True, pdf = True, pgf = False, name='series')

"""=== LOADING DATA ==="""
# check whether datafilename is given/data file exists
if len(datafilename) == 0:
    print "ERROR: No path to data file given."
    sys.exit(1)
if not os.path.exists(datafilename):
    print "ERROR: '%s' doesn't exist." % datafilename
    sys.exit(1)

# read from datafile
print "Reading data from '%s.'" % datafilename
datafile = open(datafilename, 'r')
s0, s1, s2, s3, s4  = pickle.load(datafile)
datafile.close()

"""==== DEFINITIONS ==="""

# autocorrelation function, normalized
def acf_n(x):
    N = len(x)
    out = np.zeros((N))
    for i in range(N):
        out[i] = np.mean(x*np.roll(x, i))
    out -= np.mean(x)**2
    out /= out[0]
    return out

# cross correlation via fft
def ccr(a, b):
    return np.fft.irfft(np.fft.rfft(a).conjugate()*np.fft.rfft(b))

# autocorrelation via ccr
def acf(x):
    x0 = x - np.mean(x)
    out = ccr(x0, x0)
    return out/out[0]

# estimation of autocorrelation time
def est_act(x):
    tau = 0.5
    kmax = 1
    xcor = acf(x)
    while kmax < 6*tau:
        tau += xcor[kmax]
        kmax += 1
    return tau, kmax
        
# automatic error analysis via autocorrelation analysis
def aea(x):
    N = len(x)
    tau, kmax = est_act(x)
    errtau = tau*np.sqrt(2*(2*kmax+1)/N)
    return np.mean(x), np.sqrt(np.var(x)), tau, errtau, int(12*(tau/errtau)**2)

print "Computing estimated autocorrelation time"
print aea(s0)
print aea(s1)
print aea(s2)
print aea(s3)
print aea(s4)



"""=== AUTOCORRELATE DATASETS ==="""
print "Computing normalized autocorrelation function of..."
print "... dataset 1"
acf0 = acf(s0)

print "... dataset 2"
acf1 = acf(s1)
print "... dataset 3"
acf2 = acf(s2)
print "... dataset 4"
acf3 = acf(s3)
print "... dataset 5"
acf4 = acf(s4)


"""==== PLOTTING ===="""
print "Plotting..."

# first 1000 values of the data series s0, s1, s2, s3, s4 over time
p.new(xlabel='time',ylabel='value')
p.plot(s0[0:1000], label='dataset 1')
p.plot(s1[0:1000], label='dataset 2')
p.plot(s2[0:1000], label='dataset 3')
p.plot(s3[0:1000], label='dataset 4')
p.plot(s4[0:1000], label='dataset 5')

# plott autocorrelation of s0, s1, s2, s3, s4 over k
p.new(xlabel='time',ylabel='normalized autocorrelation')
p.plot(acf0[0:1000], label='acf of dataset 1')
p.plot(acf1[0:1000], label='acf of dataset 2')
p.plot(acf2[0:1000], label='acf of dataset 3')
p.plot(acf3[0:1000], label='acf of dataset 4')
p.plot(acf4[0:1000], label='acf of dataset 5')

p.make()

print "Finished..."