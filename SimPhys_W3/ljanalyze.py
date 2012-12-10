#!/usr/bin/python2
# -*- coding:utf-8 -*-

from __future__ import division
import sys, os, pickle
import numpy as np
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
    datafilename = "data/ljsim.dat"

# check whether data file exists
if not os.path.exists(datafilename):
    print "ERROR: %s doesn't exist."
    sys.exit(1)

print "Reading data from %s." % datafilename
datafile = open(datafilename, 'r')
ts, Es, Epots, Ekins, Ts, Ps, traj = pickle.load(datafile)
datafile.close()

# drop the beginning
i = ts >=15
Es=Es[i].copy()
ts=ts[i].copy()
Epots=Epots[i].copy()
Ekins=Ekins[i].copy()
Ts=Ts[i].copy()
Ps=Ps[i].copy()
traj=traj[i].copy()

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
meanEs=np.mean(Es[:10])
print "meanEs=", meanEs
meanEpots=np.mean(Epots[:10])
print "meanEpots=", meanEpots
meanEkins=np.mean(Ekins[:10])
print "meanEkins=", meanEkins
meanTs=np.mean(Ts[:10])
print "meanTs=", meanTs
meanPs=np.mean(Ps[:10])
print "meanPs=", meanPs

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
ts10=ts[4:-5]
ts100=ts[49:-50]
#Averages
Es10=compute_running_average(Es,10)
Es100=compute_running_average(Es,100)

p.new(xlabel='time',ylabel='energy')
p.plot(ts,Es, label='Eges')
p.plot(ts10,Es10, label='running av 10')
p.plot(ts100,Es100, label='running av 100')

# Energies till time t=50
p.new(xlabel='time',ylabel='energy')
p.plot(ts[:50],Ekins[:50], label='Ekin')
p.plot(ts[:50],Es[:50], label='Eges')
p.plot(ts[:50],Epots[:50], label='Epot')

# Energies
# Averages
Ekins10=compute_running_average(Ekins,10)
Ekins100=compute_running_average(Ekins,100)
Epots10=compute_running_average(Epots,10)
Epots100=compute_running_average(Epots,100)


p.new(xlabel='time',ylabel='energy')
p.plot(ts,Ekins, label='Ekin')
p.plot(ts10,Ekins10, label='running av 10')
p.plot(ts100,Ekins100, label='running av 100')

p.plot(ts,Es, label='Eges')
p.plot(ts,Epots, label='Epot')
p.plot(ts10,Epots10,label='running av 10')
p.plot(ts100,Epots100,label='running av 100')

# Temperature
#Average
Ts10=compute_running_average(Ts, 10)
Ts100=compute_running_average(Ts, 100)

p.new(xlabel='time',ylabel='temperature')
p.plot(ts,Ts, label='T')
p.plot(ts10,Ts10, label='running av 10')
p.plot(ts100,Ts100, label='running av 100')

# Pressure
#Average
Ps10=compute_running_average(Ps, 10)
Ps100=compute_running_average(Ps, 100)
p.new(xlabel='time',ylabel='pressure')
p.plot(ts,Ps, label='P')
p.plot(ts10,Ps10, label='running av 10')
p.plot(ts100,Ps100, label='running av 100')

p.make(ncols= 3)

print "Finished."
