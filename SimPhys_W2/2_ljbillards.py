#!/usr/bin/python2
# -*- coding:utf-8 -*-

from __future__ import division
import numpy as np

# ---- this section is usefull everytime ----
import os, sys
sys.path.append('../libs')
from evaluation import Plotter
filename = os.path.splitext(os.path.basename(__file__))[0]
p = Plotter(show = True, save = False, pgf = True, name=filename, directory = '')
# -------------------------------------------

# ==== DEFINITIONS ====
EPS = 1
SIG = 1

# constants
dt = 0.01
tmax = 20.0

# running variables
t = 0.0

# particle positions
x = np.zeros((3,5))
x[:,0] = [0.0, 0.0, 0.0]
x[:,1] = [5.0, 0.3, 0.0]
x[:,2] = [8.0, 1.8, 0.0]
x[:,3] = [11.0, 0.0, -1.0]
#x[:,4] = [12.0, 9.0, 0.0]
x[:,4] = [15.4324, 9.51146, 0.0]

# particle velocities
v = np.zeros((3,5))
v[:,0] = [2.0, 0.0, 0.0]
v[:,1] = [0.0, 0.0, 0.0]
v[:,2] = [0.0, 0.0, 0.0]
v[:,3] = [0.0, 0.0, 0.0]
v[:,4] = [0.0, 0.0, 0.0]

# ==== FUNCTIONS ====
def compute_lj_potential(rij, eps = EPS, sig = SIG):
    q = sig/np.linalg.norm(rij)
    return 4*eps*(q**12-q**6)

def compute_lj_force(rij, eps = EPS, sig = SIG):
    norm = np.linalg.norm(rij)
    q = sig/norm
    return 4*eps*(12*q**11-6*q**5)*q/norm**2*rij

def compute_forces(x):
    """Compute and return the forces acting onto the particles,
    depending on the positions x."""
    global epsilon, sigma
    _, N = x.shape
    f = np.zeros_like(x)
    for i in range(1,N):
        for j in range(i):
            # distance vector
            rij = x[:,j] - x[:,i]
            fij = compute_lj_force(rij)
            f[:,i] -= fij
            f[:,j] += fij
    return f

def compute_energy(x, v):
    """Compute and return the total energy of the system with the
    particles at positions x."""
    _, N = x.shape
    E_pot = 0.0
    E_kin = 0.0
    # sum up potential energies
    for i in range(1,N):
        for j in range(i):
            # distance vector
            rij = x[:,j] - x[:,i]
            E_pot += compute_lj_potential(rij)
    # sum up kinetic energy
    for i in range(N):
        E_kin += 0.5 * np.dot(v[:,i],v[:,i])
    return E_pot + E_kin
    
def step_vv(x, v, f, dt):
    # update positions
    x += v*dt + 0.5*f * dt*dt
    # half update of the velocity
    v += 0.5*f * dt
        
    # compute new forces
    f = compute_forces(x)
    # we assume that m=1 for all particles

    # second half update of the velocity
    v += 0.5*f * dt

    return x, v, f

# ==== CALCULATION ====
f = compute_forces(x)

# variables to cumulate data
traj = []
Es = []

# number of particles
N = x.shape[1]

# open the trajectory file
vtffile = open('./plots/vmd_2_ljbillards/ljbillards.vtf', 'w')
# write the structure of the system into the file: 
# N particles ("atoms") with a radius of 0.5
vtffile.write('atom 0:%s radius 0.5\n' % (N-1))

# main loop
while t < tmax:
    x, v, f = step_vv(x, v, f, dt)
    t += dt

    traj.append(x.copy())
    Es.append(compute_energy(x, v))
    
    # write out that a new timestep starts
    vtffile.write('timestep\n')
    # write out the coordinates of the particles
    for i in range(N):
        vtffile.write("%s %s %s\n" % (x[0,i], x[1,i], x[2,i]))

vtffile.close()

traj = np.array(traj)

# ==== PLOTTING ====
# plot the trajectory
p.new(aspect='equal')
for i in range(N):
    p.plot(traj[:,0,i], traj[:,1,i], label='%s'%i)

# plot the total energy
p.new()
p.plot(Es, label='energy')

p.make()