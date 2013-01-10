#!/usr/bin/python2
# -*- coding:utf-8 -*-

from __future__ import division
import sys, os, pickle
import numpy as np
from libs.cython import set_globals, compute_distances
from libs.simlib import Plotter
from matplotlib import cm

p = Plotter(show = True, pdf = False, pgf = False, name='ljsim')

"""==== DEFINITIONS ===="""
# SYSTEM CONSTANTS
# density
density = 0.316
# number of particles per side for cubic setup
n = 10

# COMPUTED CONSTANTS
# total number of particles
N = n*n*n
# volume of the system
volume = N/density
# side length of the system
L = volume**(1./3.)

# get filename on command line
if len(sys.argv) == 2:
    print "Usage: python %s FILE" % sys.argv[0]
    datafilename = sys.argv[1]
else:
    datafilename = "data/03/ljsim.dat"

# check whether data file exists
if not os.path.exists(datafilename):
    print "ERROR: %s doesn't exist."
    sys.exit(1)

print "Reading data from %s." % datafilename
datafile = open(datafilename, 'r')
ts, Es, Epots, Ekins, Ts, Ps, traj = pickle.load(datafile)
datafile.close()

set_globals(L, N, 42, 42)

s0 = Es
s1 = Epots
s2 = Ekins
s3 = Ts
s4 = Ps

k=500

"""==== FUNCTIONS ==="""
def variance(x):
    out = (x-np.mean(x))**2
    out = np.sum(out)
    return out/len(x)/(len(x)-1)

def mean(x): return np.sum(x)/len(x)

# autocorrelation function, normalized
def acf_n(x):
    N = len(x)
    x -= mean(x)
    out = np.zeros((N))
    for i in range(N):
        out[i] = np.mean(x[i:N]*x[0:N-i])
    #out -= np.mean(x)**2
    return out/out[0]
              
# integrated autocorrelation function
def acf_int(x):
    N = len(x)
    x0 = acf(x)
    out = np.zeros((N))
    for i in range(N):
        if i == 0: out[i] = x0[0]
        else: out[i] = out[i-1]+x0[i]
    return out

# cross correlation via fft
def ccr(a, b):
    return np.fft.irfft(np.fft.rfft(a).conjugate()*np.fft.rfft(b))

# autocorrelation via ccr
def acf(x):
    x0 = x - np.mean(x)
    x0 = np.append(x0, np.zeros((len(x))))
    out = ccr(x0, x0)
    out = out [0:len(x)]
    out /= np.arange(len(x), 0, -1)
    return (out/out[0])[:len(x)]

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
    return np.mean(x), np.sqrt(variance(x)), tau, errtau, int(12*(tau/errtau)**2)

# BINNING ANALYSIS
# blocking of given time series x and block size k
def block(x, k):
    N = len(x)    
    x = x[0:N-N%k]
    N = len(x)
    B = int(N/k)
    out = np.zeros((B))
    for i in range(B):
        out[i] = np.mean(x[i*k:(i+1)*k])
    return out

# compute block variance for given time series x and block size k
def cbv(x, k):
    b = block(x, k)
    b -= np.mean(b)
    return np.sum(b*b)/(len(b)-1)

# estimated autocorrelation time from blocking/blocking tau
def est_act_b(x, k): return 0.5*k*cbv(x,k)/np.var(x)

# estimate error of mean value using block variance
def eem(x, k): return np.sqrt(cbv(x, k)/(len(x)//k))

# compute sequences of block variance for plotting
def cbv_seq(x):
    global k
    t = np.zeros((len(x)))
    for i in range(k):
        t[i] = est_act_b(x, i+1)
    return t[0:k]

# compute sequences of estimated error of mean value using block variance
def eem_seq(x):
    global k
    t = np.zeros((len(x)))
    for i in range(k):
        t[i] = eem(x, i+1)
    return t[0:k]

# WORKING, BUT VERY SLOW
# blocking for jackknifing
def block_j(x,k):
    N = len(x)    
    x = x[0:N-N%k]
    N = len(x)
    B = N//k
    out = np.zeros((B,k))
    for i in range(B):
        out[i,:] = x[i*k:(i+1)*k]
    return out

# Jackknife error
def jke(x, k):
    xb = block_j(x, k)
    #Nb = len(xb[:,0])
    Nb = np.shape(xb)[0]
    Oj = np.mean(x)
    out = 0
    for i in range(Nb):
        h = (np.append(xb[0:i,:],xb[i+1:,:]) - Oj)
        out += np.dot(h,h)
    return np.sqrt(out*(Nb-1)/Nb)

    
# compute sequences of estimated error of mean value using block variance
def jke_seq(x):
    global k
    t = np.zeros((len(x)))
    for i in range(k):
        t[i] = jke(x, i+1)
    return t[0:k]


"""=== AUTOCORRELATE DATASETS ==="""
print "Computing normalized autocorrelation function of..."
print "... total energy"
acfn0 = acf_n(s0)
print "... potential energy"
acfn1 = acf_n(s1)
print "... kinetic energy"
acfn2 = acf_n(s2)
print "... temperature"
acfn3 = acf_n(s3)
print "... pressure"
acfn4 = acf_n(s4)

print "Computing normalized autocorrelation function (via fft) of..."
print "... total energy"
acf0 = acf(s0)
print "... potential energy"
acf1 = acf(s1)
print "... kinetic energy"
acf2 = acf(s2)
print "... temperature"
acf3 = acf(s3)
print "... pressure"
acf4 = acf(s4)

print "Computing integrated autocorrelation function (via fft) of..."
print "... total energy"
acfi0 = acf_int(s0)
print "... potential energy"
acfi1 = acf_int(s1)
print "... kinetic energy"
acfi2 = acf_int(s2)
print "... temperature"
acfi3 = acf_int(s3)
print "... pressure"
acfi4 = acf_int(s4)

"""=== AUTOMATIC ERROR ANALYSIS ==="""
print "Computing automatic error analysis function of..."
print "| used dataset | mean value | error of mean value | est. autocor. time | error of autocor. time | N_eff | "
print "... total energy:", aea(s0)
print "... potential energy:", aea(s1)
print "... kinetic energy:", aea(s2)
print "... temperature:", aea(s3)
print "... pressure:", aea(s4)

"""=== BINNING ANALYSIS ==="""
print "Computing blocking taus of..."
print "... total energy"
bts0 = cbv_seq(s0)
print "... potential energy"
bts1 = cbv_seq(s1)
print "... kinetic energy"
bts2 = cbv_seq(s2)
print "... temperature"
bts3 = cbv_seq(s3)
print "... pressure"
bts4 = cbv_seq(s4)

print "Computing blocking estimated error of mean value of..."
print "... total energy"
ems0 = eem_seq(s0)
print "... potential energy"
ems1 = eem_seq(s1)
print "... kinetic energy"
ems2 = eem_seq(s2)
print "... temperature"
ems3 = eem_seq(s3)
print "... pressure"
ems4 = eem_seq(s4)

"""=== JACKKNIFE ANALYSIS ==="""
print "Computing jackknife error of mean value of..."
print "... total energy"
jks0 = jke_seq(s0)
print "... potential energy"
jks1 = jke_seq(s1)
print "... kinetic energy"
jks2 = jke_seq(s2)
print "... temperature"
jks3 = jke_seq(s3)
print "... pressure"
jks4 = jke_seq(s4)


"""==== PLOTTING ===="""
print "Plotting..."

# first 1000 values of the data series s0, s1, s2, s3, s4 over time
p.new(xlabel='time',ylabel='value')
p.plot(s0, label='total energy')
p.plot(s1, label='potential energy')
p.plot(s2, label='kinetic energy')
p.plot(s3, label='temperature')
p.plot(s4, label='pressure')

s0 = Es
s1 = Epots
s2 = Ekins
s3 = Ts
s4 = Ps



# plot autocorrelation of s0, s1, s2, s3, s4 over k
p.new(xlabel='k',ylabel='normalized autocorrelation')
p.plot(acfn0, label='acf of total energy')
p.plot(acfn1, label='acf of potential energy')
p.plot(acfn2, label='acf of kinetic energy')
p.plot(acfn3, label='acf of temperature')
p.plot(acfn4, label='acf of pressure')

# plot autocorrelation via fft  of s0, s1, s2, s3, s4 over k
p.new(xlabel='k',ylabel='normalized autocorrelation via fft')
p.plot(acf0, label='acf of total energy')
p.plot(acf1, label='acf of potential energy')
p.plot(acf2, label='acf of kinetic energy')
p.plot(acf3, label='acf of temperature')
p.plot(acf4, label='acf of pressure')

# plot autocorrelation via fft  of s0, s1, s2, s3, s4 over k, ZOOM to relevant interval
p.new(xlabel='k',ylabel='normalized autocorrelation via fft, ZOOM')
p.plot(acf0[0:100], label='acf of total energy')
p.plot(acf1[0:100], label='acf of potential energy')
p.plot(acf2[0:100], label='acf of kinetic energy')
p.plot(acf3[0:100], label='acf of temperature')
p.plot(acf4[0:100], label='acf of pressure')

# plot integrated autocorrelation via fft  of s0, s1, s2, s3, s4 over k
p.new(xlabel='k',ylabel='integrated autocorrelation function via fft')
p.plot(acfi0, label='acf of total energy')
p.plot(acfi1, label='acf of potential energy')
p.plot(acfi2, label='acf of kinetic energy')
p.plot(acfi3, label='acf of temperature')
p.plot(acfi4, label='acf of pressure')

# plot integrated autocorrelation via fft  of s0, s1, s2, s3, s4 over k, ZOOM to relevant interval
p.new(xlabel='k',ylabel='integrated autocorrelation function via fft, ZOOM')
p.plot(acfi0[0:100], label='acf of total energy')
p.plot(acfi1[0:100], label='acf of potential energy')
p.plot(acfi2[0:100], label='acf of kinetic energy')
p.plot(acfi3[0:100], label='acf of temperature')
p.plot(acfi4[0:100], label='acf of pressure')

# plot blocking tau over block size k
p.new(xlabel='block size k',ylabel='blocking tau')
p.plot(bts0, label='bt of total energy')
p.plot(bts1, label='bt of potential energy')
p.plot(bts2, label='bt of kinetic energy')
p.plot(bts3, label='bt of temperature')
p.plot(bts4, label='bt of pressure')

# plot blocking eem over block size k
p.new(xlabel='block size k',ylabel='est. err. of mean value')
p.plot(ems0, label='eem of total energy')
p.plot(ems1, label='eem of potential energy')
p.plot(ems2, label='eem of kinetic energy')
p.plot(ems3, label='eem of temperature')
p.plot(ems4, label='eem of pressure')

# plot jke over block size k
p.new(xlabel='block size k',ylabel='jackknife error')
p.plot(jks0, label='jke of total energy')
p.plot(jks1, label='jke of potential energy')
p.plot(jks2, label='jke of kinetic energy')
p.plot(jks3, label='jke of temperature')
p.plot(jks4, label='jke of pressure')

p.make(ncols=3,savewindow=True)

print "Finished..."