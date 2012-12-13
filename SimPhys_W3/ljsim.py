#!/usr/bin/python2
# -*- coding:utf-8 -*-

from __future__ import division
import os, pickle
import numpy as np
from libs.cython import set_globals, compute_forces, compute_distances, compute_energy, compute_pressure, rebuild_neighbor_lists
from libs.simlib import Plotter
from matplotlib import cm

p = Plotter(show = True, pdf = False, pgf = False, name='ljsim')

"""==== DEFINITIONS ===="""
# SYSTEM CONSTANTS
# density
density = 0.316
# timestep
dt = 0.01
# max length of each run 800
tadd = 20.0
# max length of all runs
tges = 1000.0
# number of particles per side for cubic setup
n = 10
# Desired temperature {0.3, 1.0, 2.0}
Tdes= 10.0

NEWRUN = True
RANDOMPOSITION = True
FORCECAPPING = True
VELOCITYRESCALING = True

# SIMULATION CONSTANTS
# skin size
skin = 0.4
# number of steps to do before the next measurement
measurement_stride = 100
# cutoff length
rcut = 2.5
# potential shift
shift = -0.016316891136
# force limit
fmax = 10
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
    global rcut, skin, fmax

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
    f = compute_forces(x,fmax)
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
    
def compute_temperature(Ekin, N = N):    
    return 2/3*Ekin/N # T_red_units = kB*T/eps, eps = 1

"""==== SIMULATION ===="""
# ==== INITIALIZATION ====
if NEWRUN:
    if os.path.exists(vtffilename): os.remove(vtffilename)
    if os.path.exists(datafilename): os.remove(datafilename)
    if os.path.exists(statefilename): os.remove(statefilename)

# SET UP SYSTEM OR LOAD IT
if os.path.exists(statefilename):
    step, t, x, v, fmax = np.load(statefilename)
    
    print "Old data was found. Restarting simulation at t=%s, step=%s with density=%s, L=%s, N=%s." % (t,step,density, L, N)
else:
    step = 0
    t = 0.0
    if RANDOMPOSITION:
        x = L*np.random.random((3,N))
    else:
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
    
    # random particle velocities -> new
    v = 0.1*(2.0*np.random.random((3,N))-1.0)

    print "No old data was found. Starting simulation with density=%s, L=%s, N=%s." %(density, L, N)

if not FORCECAPPING: fmax = 0

# calculate number of steps
tadd = min(tges-t, tadd)
steps = int(tadd//(dt*measurement_stride))

# variables to cumulate data
traj = np.empty((steps,3,N))
ts = np.empty(steps)
Es = np.empty(steps)
Epots = np.empty(steps)
Ekins = np.empty(steps)
Ts = np.empty(steps)
Ps = np.empty(steps)

# ==== CALCULATION ====
print "Simulating until tmax=%s..." % (t + tadd)

set_globals(L, N, rcut, shift)
rebuild_neighbor_lists(x, rcut+skin)
xup = x.copy()
f = compute_forces(x,fmax)

# calculate or load the data from the time before the current run will start
if os.path.exists(datafilename):
    ts_old, Es_old, Epots_old, Ekins_old, Ts_old, Ps_old, traj_old = np.load(datafilename)
else:
    a,b,c = compute_energy(x, v)
    d = compute_temperature(c)
    e = compute_pressure(x,v)
    ts_old = np.array([t])
    Es_old = np.array([a])
    Epots_old = np.array([b])
    Ekins_old = np.array([c])
    Ts_old = np.array([d])
    Ps_old = np.array([e])
    traj_old = np.array([x])

# main loop
for n in range(steps):
    for _i in range(measurement_stride):
        x, v, f, xup = step_vv(x, v, f, dt, xup)
        t += dt
        step += 1
    E, Epot, Ekin = compute_energy(x,v)
    T = compute_temperature(Ekin)
    P = compute_pressure(x,v)
    print "t=%s, E=%s, Epot=%s, Ekin=%s, T=%s, P=%s" % (t, E, Epot, Ekin,T,P)
    
    #Velocity rescaling
    if VELOCITYRESCALING: v*=np.sqrt(Tdes/T)
    
    # store data
    ts[n] = t
    Es[n] = E
    Epots[n] = Epot
    Ekins[n] = Ekin
    Ts[n] = T
    Ps[n] = P
    traj[n] = x
    fmax *= 1.1
    
print "Finished simulation."

"""==== SAVING ===="""
# write out vtf
print "Writing vtf data to %s ..." % vtffilename
write_vtf(traj)

# write out state data
print "Writing state data to %s ..." % statefilename
statefile = open(statefilename, 'w')
pickle.dump([step, t, x, v, fmax], statefile)
statefile.close()

# write out simulation data
print "Writing simulation data to %s ..." % datafilename
ts = np.append(ts_old,ts,axis = 0)
Es = np.append(Es_old,Es,axis = 0)
Epots = np.append(Epots_old,Epots,axis = 0)
Ekins = np.append(Ekins_old,Ekins,axis = 0)
Ts = np.append(Ts_old,Ts,axis = 0)
Ps = np.append(Ps_old,Ps,axis = 0)
traj = np.append(traj_old,traj,axis = 0)
datafile = open(datafilename, 'w')
pickle.dump([ts, Es, Epots, Ekins, Ts, Ps, traj], datafile)
datafile.close()

"""==== PLOTTING ===="""
print "Plotting..."

# Trajectories
p.new(aspect='equal',xlabel='x-coordinate',ylabel='y-coordinate')
traj -= np.floor(traj/L)*L

p.plot([0,L,L,0,0],[0,0,L,L,0],'b-', lw=2)
p.plot(traj[-1,0,:],traj[-1,1,:],'wo', alpha=0.1 ,ms=7, mew = 2)
p.plot(traj[0,0,:],traj[0,1,:],'+', c=[0.8,0.8,0.8], alpha=0.5)

i = range(traj.shape[2])
np.random.shuffle(i)
tpart = np.array(traj[:,:,i[:3]])
for n in range(1,tpart.shape[0]):
    i = (tpart[n-1,0,:] - tpart[n,0,:])**2+(tpart[n-1,1,:] - tpart[n,1,:])**2 > L*L*0.5
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

# RDF
datafile = open(datafilename, 'r')
ts, Es, Epots, Ekins, Ts, Ps, traj2 = pickle.load(datafile)
datafile.close()
p.new(xlabel='distance',ylabel='probability')
p.hist(compute_distances(traj2[-1]), bins=100, range=(0.8,5),normed=True,  log=False, label='RDF')

p.make(ncols= 3)

print "Finished."
