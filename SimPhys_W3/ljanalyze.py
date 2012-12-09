#!/usr/bin/python2
# -*- coding:utf-8 -*-

from __future__ import division
import sys, os, pickle
import numpy as np
from libs.simlib import Plotter
from matplotlib import cm

p = Plotter(show = True, pdf = False, pgf = False, name='ljfluid')

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
    datafilename = "data/until1000/ljsim.dat"

# check whether data file exists
if not os.path.exists(datafilename):
    print "ERROR: %s doesn't exist."
    sys.exit(1)

print "Reading data from %s." % datafilename
datafile = open(datafilename, 'r')
ts, Es, Epots, Ekins, Ts, Ps, traj = pickle.load(datafile)
datafile.close()

"""==== PLOTTING ===="""
print "Plotting..."

# Trajectories
p.new(aspect='equal',xlabel='x-coordinate',ylabel='y-coordinate')
traj -= np.floor(traj/L)*L

p.plot([0,L,L,0,0],[0,0,L,L,0],'b-', lw=2)
p.plot(traj[-1,0,:],traj[-1,1,:],'wo', alpha=0.1 ,ms=7, mew = 2)
p.plot(traj[0,0,:],traj[0,1,:],'+', c=[0.8,0.8,0.8], alpha=0.1)

i = range(traj.shape[2])
np.random.shuffle(i)
tpart = np.array(traj[:,:,i[:3]])
for n in range(1,tpart.shape[0]):
    i = (tpart[n-1,0,:] - tpart[n,0,:])**2+(tpart[n-1,1,:] - tpart[n,1,:])**2 > 50
    tpart[n,:,i] = [None,None,None]
nmax = tpart.shape[2]
cm = cm.get_cmap('Dark2')
colors=[cm(1.*i/nmax) for i in range(nmax)]
for n in range(nmax):
    p.plot(tpart[:,0,n],tpart[:,1,n],'-', c = colors[n], alpha=0.8)
    p.plot(tpart[-1,0,n],tpart[-1,1,n],'o', c = colors[n], alpha=0.8 ,ms=7, mew = 2)

# Total energy
p.new(xlabel='time',ylabel='energy')
p.plot(ts,Es, label='Eges')

# Energies
p.new(xlabel='time',ylabel='energy')
p.plot(ts,Ekins, label='Ekin')
p.plot(ts,Es, label='Eges')
p.plot(ts,Epots, label='Epot')

# Temperature
p.new(xlabel='time',ylabel='temperature')
p.plot(ts,Ts, label='T')

# Pressure
p.new(xlabel='time',ylabel='pressure')
p.plot(ts,Ps, label='P')

p.make(ncols= 3)

print "Finished."