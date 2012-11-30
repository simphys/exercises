#!/usr/bin/python2
# -*- coding:utf-8 -*-

from __future__ import division
from numpy import *
from libs.cython2.lj import *
from libs.simlib import Plotter
import sys

p = Plotter(show = True, pdf = False, pgf = False, name='7_ljfluid_lists')


# CONSTANTS
# density
density = 0.7
# number of particles per side
if len(sys.argv) == 2 and sys.argv[1].isdigit(): n = int(sys.argv[1])
else: n = 5
# timestep
dt = 0.01
# length of run
tmax = 50.0
# cutoff length
rcut = 2.5
# potential shift
shift = -0.016316891136
# skin size
skin = 0.3

# COMPUTED CONSTANTS
# total number of particles
N = n*n*n
# volume of the system
volume = N/density
# side length of the system
L = volume**(1./3.)

# RUNNING VARIABLES
t = 0.0

# particle positions on cubic lattice
x = empty((3,N))
l = L/n
count = 0
for i in range(n):
    for j in range(n):
        for k in range(n):
            x[:,count] = [i*l, j*l, k*l]
            count += 1

# random particle velocities
random.seed(17)
v = 2.0*random.random((3,N))-1.0

# variables to cumulate data
ts = []
Es = []

# open the trajectory file
vtffile = open('./plots/vmd_4_ljfluid/ljfluid_lists.vtf', 'w')
# write the structure of the system into the file: 
# N particles ("atoms") with a radius of 0.5
vtffile.write('atom 0:%s radius 0.5\n' % (N-1))
vtffile.write('pbc %s %s %s\n' % (L, L, L))

# write out that a new timestep starts
vtffile.write('timestep\n')
# write out the coordinates of the particles
for i in range(N):
    vtffile.write("%s %s %s\n" % (x[0,i], x[1,i], x[2,i]))

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

# main loop
set_globals(L, N, rcut, shift)
rebuild_neighbor_lists(x, rcut+skin)
xup = x.copy()
f = compute_forces(x)

traj = np.array([x[0:2,:] - np.floor(x[0:2,:]/L)*L])
while t < tmax:
    x, v, f, xup = step_vv(x, v, f, dt, xup)
    t += dt

    E = compute_energy(x, v)
    
    print "t=%s, E=%s" % (t, E)

    ts.append(t)
    Es.append(E)

    # write out that a new timestep starts
    vtffile.write('timestep\n')
    # write out the coordinates of the particles
    for i in range(N):
        vtffile.write("%s %s %s\n" % (x[0,i], x[1,i], x[2,i]))
    
    # for plotting
    traj = np.append(traj,[x[0:2,:] - np.floor(x[0:2,:]/L)*L],axis=0)

vtffile.close()

p.new()
for n in range(1,traj.shape[0]):
    i = (traj[n-1,0,:] - traj[n,0,:])**2+(traj[n-1,1,:] - traj[n,1,:])**2 > 1
    traj[n,:,i] = [None,None]
p.plot(traj[:,0,:],traj[:,1,:],'-', lw = 2, alpha=0.3)
p.plot(traj[0,0,:],traj[0,1,:],'ko', alpha=0.5 ,ms=10, mew = 4, label = 'start position')
p.plot(traj[-1,0,:],traj[-1,1,:],'wo', alpha=0.5 ,ms=10, mew = 4, label = 'end position')

p.new()
p.plot(ts, Es)
p.make(ncols=1)
