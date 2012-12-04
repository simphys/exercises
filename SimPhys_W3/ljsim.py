#!/usr/bin/python2
# -*- coding:utf-8 -*-

from __future__ import division
import os, pickle
import numpy as np
from libs.cython import set_globals, compute_forces, compute_energy, rebuild_neighbor_lists
from libs.simlib import Plotter
from matplotlib import cm

p = Plotter(show = True, pdf = False, pgf = False, name='ljfluid')

"""==== DEFINITIONS ===="""
# SYSTEM CONSTANTS
# density
density = 0.316
# timestep
dt = 0.01
# max length of each run
tadd = 50.0
# max length of all runs
tges = 1000.0
# number of particles per side for cubic setup
n = 10

# SIMULATION CONSTANTS
# skin size
skin = 0.4
# number of steps to do before the next measurement
measurement_stride = 100
# cutoff length
rcut = 2.5
# potential shift
shift = -0.016316891136
# VTF filename 
vtffilename = "ljsim.vtf"
# DATA filename 
datafilename = "ljsim.dat"

# COMPUTED CONSTANTS
# total number of particles
N = n*n*n
# volume of the system
volume = N/density
# side length of the system
L = volume**(1./3.)

# SEED
np.random.seed(42)

"""==== FUNCTIONS ===="""
def step_vv(x, v, f, dt, xup):
    global rcut, skin

    # update positions
    x += v*dt + 0.5*f * dt*dt

    # compute maximal position update
    # vectorial
    dx = x - xup
    # square
    dx *= dx
    # sum up 
    dx = dx.sum(axis=0)
    # test whether the neighbor list needs to be rebuilt
    if max(dx) > (0.5*skin)**2:
        rebuild_neighbor_lists(x, rcut+skin)
        xup = x.copy()
    
    # half update of the velocity
    v += 0.5*f * dt
        
    # compute new forces
    f = compute_forces(x)
    # we assume that m=1 for all particles

    # second half update of the velocity
    v += 0.5*f * dt

    return x, v, f, xup

def write_vtf(traj_new, vtffilename = vtffilename):    
    # check whether vtf file already exists
    if os.path.exists(vtffilename):
        print "Opening %s to append new timesteps." % vtffilename
        vtffile = open(vtffilename, 'a')
    else:
        print "Creating %s..." % vtffilename
        # create a new file and write the structure
        vtffile = open(vtffilename, 'w')
    
        # write the structure of the system into the file: 
        # N particles ("atoms") with a radius of 0.5
        vtffile.write('atom 0:%s radius 0.5\n' % (N-1))
        vtffile.write('pbc %s %s %s\n' % (L, L, L))
    
    # write vtf file
    for x in traj_new:
        # write out that a new timestep starts
        vtffile.write('timestep\n')
        # write out the coordinates of the particles
        for i in range(N):
            vtffile.write("%s %s %s\n" % (x[0,i], x[1,i], x[2,i]))
        
    # close vtf file
    print "Closing %s." % vtffilename
    vtffile.close()

"""==== SIMULATION ===="""
# ==== INITIALIZATION ====
# SET UP SYSTEM OR LOAD IT
# check whether data file already exists
if os.path.exists(datafilename):
    print "Reading data from %s." % datafilename
    step, t, x, v, ts, Es, traj = np.load(datafilename)
    print "Restarting simulation at t=%s..." % t

else:
    print "Starting simulation..."
    t = 0.0
    step = 0
    # particle positions on cubic lattice
    x = np.empty((3,N))
    count = 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                x[:,count] = [i, j, k]
                count += 1
    x += 0.5
    x *= L/n
    
    # random particle velocities
    v = 0.1*(2.0*np.random.random((3,N))-1.0)
    
    # start values
    ts = np.array([])
    Es = np.array([])
    traj = np.array([x])

# calculate number of steps
tadd = min(tges-t, tadd)
steps = int(tadd//(dt*measurement_stride))

# variables to cumulate data
traj_new = np.empty((steps,3,N))
ts_new = np.empty(steps)
Es_new = np.empty(steps)

print "density=%s, L=%s, N=%s" % (density, L, N)

# ==== CALCULATION ====
# main loop
set_globals(L, N, rcut, shift)
rebuild_neighbor_lists(x, rcut+skin)
xup = x.copy()
f = compute_forces(x)

print "Simulating until tmax=%s..." % (t + tadd)

for n in range(steps):
    for _i in range(measurement_stride):
        x, v, f, xup = step_vv(x, v, f, dt, xup)
        t += dt
        step += 1
        
    E = compute_energy(x, v)
    print "t=%s, E=%s" % (t, E)
    
    # store data
    ts_new[n] = t
    Es_new[n] = E
    traj_new[n] = x

ts = np.append(ts,ts_new,axis = 0)
Es = np.append(Es,Es_new,axis = 0)
traj = np.append(traj,traj_new,axis = 0)

"""==== SAVING ===="""
# write out vtf
write_vtf(traj_new)

# write out simulation data
print "Writing simulation data to %s." % datafilename
datafile = open(datafilename, 'w')
pickle.dump([(n-1)*measurement_stride, t, x, v, ts, Es, traj], datafile)
datafile.close()

print "Finished simulation."

"""==== PLOTTING ===="""
print "Plotting..."

# Trajectories
p.new(aspect='equal')
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

# Energy
p.new()
p.plot(ts,Es)
p.make(ncols= 2)

print "Finished."