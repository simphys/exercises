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
    datafilename = "data/10/ljsim.dat"

# check whether data file exists
if not os.path.exists(datafilename):
    print "ERROR: %s doesn't exist."
    sys.exit(1)

print "Reading data from %s." % datafilename
datafile = open(datafilename, 'r')
ts, Es, Epots, Ekins, Ts, Ps, traj = pickle.load(datafile)
datafile.close()

set_globals(L, N, 42, 42)

"""==== Functions==== """
def compute_running_average(O,M):
    length=np.shape(O)[0]
    av=[]
    for t in range(M-1,length):
        summe = 0.0
        for i in range(1,M+1):
            summe+=O[t-M+i]
        av.append(summe/M)
    return np.array(av)

"""==== PLOTTING ===="""
# mean values
i = ts>700
meanEs=np.mean(Es[i])
print "meanEs=", meanEs
meanEpots=np.mean(Epots[i])
print "meanEpots=", meanEpots
meanEkins=np.mean(Ekins[i])
print "meanEkins=", meanEkins
meanTs=np.mean(Ts[i])
print "meanTs=", meanTs
meanPs=np.mean(Ps[i])
print "meanPs=", meanPs

print "Plotting..."

# Trajectories
p.new(title='Trajectories', aspect='equal',xlabel='x-coordinate',ylabel='y-coordinate')
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

# Energies
p.new(title='Energies (from t=0 to t=10)',xlabel='time',ylabel='energy')
i = ts < 10
p.plot(ts[i],Ekins[i], label='Ekin')
p.plot(ts[i ],Es[i], label='Eges')
p.plot(ts[i ],Epots[i], label='Epot')

# Energies
p.new(title='Energies (from t=10 to t=100)',xlabel='time',ylabel='energy')
i = np.all([ts < 100, ts > 10], axis=0)
p.plot(ts[i],Ekins[i], label='Ekin')
p.plot(ts[i ],Es[i], label='Eges')
p.plot(ts[i ],Epots[i], label='Epot')


# drop the beginning
i = ts >100
if sum(i) > 20:
    Es=Es[i].copy()
    ts=ts[i].copy()
    Epots=Epots[i].copy()
    Ekins=Ekins[i].copy()
    Ts=Ts[i].copy()
    Ps=Ps[i].copy()
    traj=traj[i].copy()

# Energies
p.new(title='Energies (from t=100)',xlabel='time',ylabel='energy')
p.plot(ts,Ekins, label='Ekin')
p.plot(ts,Es, label='Eges')
p.plot(ts,Epots, label='Epot')

# Energies
# Averages
ts10=ts[4:-5]
ts100=ts[49:-50]
Es10=compute_running_average(Es,10)
Es100=compute_running_average(Es,100)
Ekins10=compute_running_average(Ekins,10)
Ekins100=compute_running_average(Ekins,100)
Epots10=compute_running_average(Epots,10)
Epots100=compute_running_average(Epots,100)

p.new(title='Total energy (from t=100)',xlabel='time',ylabel='energy')
p.plot(ts,Es, label='Eges')
p.plot(ts10,Es10, label='running av 10')
p.plot(ts100,Es100, label='running av 100')

p.new(title='Kinetic energy (from t=100)',xlabel='time',ylabel='energy')
p.plot(ts,Ekins, label='Ekin')
p.plot(ts10,Ekins10, label='running av 10')
p.plot(ts100,Ekins100, label='running av 100')

p.new(title='Potential energy (from t=100)',xlabel='time',ylabel='energy')
p.plot(ts,Epots, label='Epot')
p.plot(ts10,Epots10,label='running av 10')
p.plot(ts100,Epots100,label='running av 100')

# Temperature
#Average
Ts10=compute_running_average(Ts, 10)
Ts100=compute_running_average(Ts, 100)

p.new(title='Temperature (from t=100)',xlabel='time',ylabel='temperature')
p.plot(ts,Ts, label='T')
p.plot(ts10,Ts10, label='running av 10')
p.plot(ts100,Ts100, label='running av 100')

# Pressure
#Average
Ps10=compute_running_average(Ps, 10)
Ps100=compute_running_average(Ps, 100)
p.new(title='Pressure (from t=100)', xlabel='time',ylabel='pressure')
p.plot(ts,Ps, label='P')
p.plot(ts10,Ps10, label='running av 10')
p.plot(ts100,Ps100, label='running av 100')

# RDF
#Average
l = min(len(traj),100)
m = N*(N-1)/2
dist = np.empty(l*m)
for i in range(1,l+1):
    dist[(i-1)*m:i*m] = compute_distances(traj[-i])  
p.new(title='RDF (last 100 trajectories)', xlabel='distance',ylabel='probability')
p.hist(dist, bins=100, range=(0.8,5), normed=True,  log=False, label='RDF')

p.make(ncols= 2, savewindow=True)

print "Finished."