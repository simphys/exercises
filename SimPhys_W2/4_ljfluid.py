#!/usr/bin/python2
# -*- coding:utf-8 -*-

from __future__ import division
import numpy as np
import sys
from libs.simlib import Plotter

p = Plotter(show = True, pdf = False, pgf = False, name='4_ljfluid')

# ==== DEFINITIONS ====
# CONSTANTS
EPS = 1
SIG = 1
# density
density = 0.7
# number of particles per side
if len(sys.argv) == 2 and sys.argv[1].isdigit(): n = int(sys.argv[1])
else: n = 3
# timestep
dt = 0.01
# length of run
tmax = 1.0
# cutoff
rcut = 2.5*SIG
# potential shift
shift = -0.016316891136
# total number of particles
N = n*n*n

# particle positions on a cubic lattice
x = np.zeros((3,N))

# compute system size
L = n/density**(1/3)

# - set up n*n*n particles on a cubic lattice
'''
# This can be done shorter, take a look below
positions = np.arange(0.5, n + 0.5)/density**(1/3)
for a in range(n):
    for b in range(n):
        for c in range(n):
            x[0,(a*n*n)+(b*n)+c] = positions[a]
            x[1,(a*n*n)+(b*n)+c] = positions[b]
            x[2,(a*n*n)+(b*n)+c] = positions[c]
'''
i = 0
for a in range(n):
    for b in range(n):
        for c in range(n):
            x[:,i] = [a+.5, b+.5, c+.5]
            i += 1
x = x/density**(1/3)
'''
p.new()
p.plot(np.array(x[0]),np.array(x[2]),'o')
p.plot(np.array(x[0]),np.array(x[1]),'o')
'''

# random particle velocities
np.random.seed(17)
v = 2.0*np.random.random((3,N))-1.0

# RUNNING VARIABLES
t = 0.0

# variables to cumulate data
ts = []
Es = []
traj = []

# ==== FUNCTIONS ====
VCOFF = 4*EPS*((SIG/rcut)**12-(SIG/rcut)**6)
def compute_lj_potential(rij, eps = EPS, sig = SIG, rcoff = rcut, vcoff = VCOFF):
    norm = np.linalg.norm(rij)
    if norm > rcoff: return 0
    q = sig/norm
    return 4*eps*(q**12-q**6) - vcoff

def compute_lj_force(rij, eps = EPS, sig = SIG, rcoff = rcut):
    norm = np.linalg.norm(rij)
    if norm > rcoff: return np.zeros(3)
    q = sig/norm
    return 4*eps*(12*q**11-6*q**5)*q/norm**2*rij

def minimum_image(ri, rj):
    global L
    # compute distance
    rij = rj-ri
    # wrap the distance into [-0.5, 0.5]
    rij -= np.rint(rij/L)*L
    return rij

def compute_forces(x, L):
    """Compute and return the forces acting onto the particles,
    depending on the positions x."""
    _, N = x.shape
    f = np.zeros_like(x)
    #rijs = np.empty_like(x)
    for i in range(1,N):
        for j in range(i):
            # distance vector
            rij = minimum_image(x[:,i], x[:,j])
            fij = compute_lj_force(rij)
            f[:,i] -= fij
            f[:,j] += fij
    return f

def compute_energy(x, v, L):
    """Compute and return the total energy of the system with the
    particles at positions x."""
    _, N = x.shape
    E_pot = 0.0
    E_kin = 0.0
    # sum up potential energies
    for i in range(1,N):
        for j in range(i):
            # distance vector
            rij = minimum_image(x[:,i], x[:,j])
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
    f = compute_forces(x, L)
    # we assume that m=1 for all particles

    # second half update of the velocity
    v += 0.5*f * dt

    return x, v, f

# open the trajectory file
vtffile = open('./plots/vmd_4_ljfluid/ljfluid.vtf', 'w')
# write the structure of the system into the file: 
# N particles ("atoms") with a radius of 0.5
vtffile.write('atom 0:%s radius 0.5\n' % (N-1))
vtffile.write('pbc %s %s %s\n' % (L, L, L))

# main loop
f = compute_forces(x, L)
while t < tmax:
    x, v, f = step_vv(x, v, f, dt)
    t += dt

    E = compute_energy(x, v, L)
    print "t=%s, E=%s" % (t, E)

    ts.append(t)
    Es.append(E)
    
    # write out that a new timestep starts
    vtffile.write('timestep\n')
    # write out the coordinates of the particles
    for i in range(N):
        vtffile.write("%s %s %s\n" % (x[0,i], x[1,i], x[2,i]))

vtffile.close()

p.new()
p.plot(ts, Es)
p.make(ncols=1)