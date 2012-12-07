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
tadd = 10.0
# max length of all runs
tges = 100.0
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
vtffilename = "data/ljsim.vtf"
# DATA filename 
datafilename = "data/ljsim.dat"
# STATE filename 
statefilename = "data/ljsim.state"

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
        vtffile = open(vtffilename, 'a')
    else:
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
    vtffile.close()

"""==== SIMULATION ===="""
# ==== INITIALIZATION ====
# SET UP SYSTEM OR LOAD IT
if os.path.exists(statefilename):
    step, t, x, v = np.load(statefilename)
    
    print "Old data was found. Restarting simulation at t=%s, step=%s with density=%s, L=%s, N=%s." % (t,step,density, L, N)
else:
    step = 0
    t = 0.0
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
    v = 1*(2.0*np.random.random((3,N))-1.0)

    print "No old data was found. Starting simulation with density=%s, L=%s, N=%s." %(density, L, N)

# calculate number of steps
tadd = min(tges-t, tadd)
steps = int(tadd//(dt*measurement_stride))

# variables to cumulate data
traj = np.empty((steps,3,N))
ts = np.empty(steps)
Es = np.empty(steps)
Epots = np.empty(steps)
Ekins = np.empty(steps)

# ==== CALCULATION ====
print "Simulating until tmax=%s..." % (t + tadd)

set_globals(L, N, rcut, shift)
rebuild_neighbor_lists(x, rcut+skin)
xup = x.copy()
f = compute_forces(x)

# calculate or load the data from the time before the current run will start
if os.path.exists(datafilename):
    ts_old, Es_old, Epots_old, Ekins_old, traj_old = np.load(datafilename)
else:
    a,b,c= compute_energy(x, v)
    ts_old = np.array([t])
    Es_old = np.array([a])
    Epots_old = np.array([b])
    Ekins_old = np.array([c])
    traj_old = np.array([x])

# main loop
for n in range(steps):
    for _i in range(measurement_stride):
        x, v, f, xup = step_vv(x, v, f, dt, xup)
        t += dt
        step += 1 
    E, Epot, Ekin = compute_energy(x, v)
    print "t=%s, E=%s, Epot=%s, Ekin=%s" % (t, E, Epot, Ekin)
    
    # store data
    ts[n] = t
    Es[n] = E
    Epots[n] = Epot
    Ekins[n] = Ekin
    traj[n] = x
    
print "Finished simulation."

"""==== SAVING ===="""
# write out vtf
print "Writing vtf data to %s ..." % vtffilename
write_vtf(traj)

# write out state data
print "Writing state data to %s ..." % statefilename
statefile = open(statefilename, 'w')
pickle.dump([step, t, x, v], statefile)
statefile.close()

# write out simulation data
print "Writing simulation data to %s ..." % datafilename
ts = np.append(ts_old,ts,axis = 0)
Es = np.append(Es_old,Es,axis = 0)
Epots = np.append(Epots_old,Epots,axis = 0)
Ekins = np.append(Ekins_old,Ekins,axis = 0)
traj = np.append(traj_old,traj,axis = 0)
datafile = open(datafilename, 'w')
pickle.dump([ts, Es, Epots, Ekins, traj], datafile)
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
p.plot(ts,Es)

# Energies
p.new(xlabel='time',ylabel='energy')
p.plot(ts,Es)
p.plot(ts,Ekins)
p.plot(ts,Epots)

p.make(ncols= 2)

print "Finished."