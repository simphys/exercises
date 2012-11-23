#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np

# ---- this section is usefull everytime ----
import os, sys
sys.path.append('../libs')
from evaluation import Plotter
filename = os.path.splitext(os.path.basename(__file__))[0]
p = Plotter(show = True, save = False, pgf = False, name=filename, directory = '')
# -------------------------------------------

name, x_init, v_init, m, g = np.load('solar_system.npy')

#print name
#print x_init
#print v_init
#print m
#print g

def compute_forces(x, m):
    global g
    _l, c = x.shape
    F = np.zeros((c, c, 2))
    for i in range(c):
        for j in range (c):
            if j > i:
                F[i,j,:] = g*m[i]*m[j]*(x[:,j] - x[:,i]) / np.linalg.norm(x[:,j] - x[:,i])**3
    F[:,:,0] -= np.transpose(F[:,:,0])
    F[:,:,1] -= np.transpose(F[:,:,1])
    Fout = np.zeros((2, c))
    for k in range(c):
        Fout[0,k] = F[k,:,0].sum()
        Fout[1,k] = F[k,:,1].sum()
    return Fout[:,:]

def step_euler(x, v, dt=.1):
    global m
    F = compute_forces(x, m)
    x += v*dt
    for i in range(len(v[0,:])):
        v[:,i] += F[:,i] / m[i] * dt   
    return x, v

def step_eulersym(x, v, dt=.1):
    global m
    F = compute_forces(x, m)
    for i in range(len(v[0,:])):
        v[:,i] += F[:,i] / m[i] * dt
    x += v*dt  
    return x, v

def step_vv(x, v, a, dt=.1):
    global m
    x += v*dt + a * 0.5 * dt*dt
    at = compute_forces(x, m)
    for i in range(len(v[0,:])):
        at[:,i] /= m[i]
    v += (a + at)/2.* dt
    return x, v, at
    
def simsys(f, dt=0.001, tend=1, xin=x_init, vin=v_init):
    t=0.
    x = np.copy(x_init)
    v = np.copy(v_init)
    traj = [np.copy(x_init)]
    while t <= tend:
        x, v = f(x, v, dt)
        traj.append(x.copy())
        t += dt
    return np.array(traj)

def simsysvv(dt=0.001, tend=1, f = step_vv, xin=x_init, vin=v_init):
    # initializie variables
    global m
    t=0.  
    x = np.copy(x_init)
    v = np.copy(v_init)
    # calculate a_0
    a = compute_forces(x, m)
    for i in range(len(v[0,:])):
        a[:,i] /= m[i]       
    traj = [np.copy(x_init)]
    while t <= tend:
        x, v, a = f(x, v, a, dt)
        traj.append(x.copy())
        t += dt
    return np.array(traj)

def printsys(traj, name, title_in):
    p.new(title=title_in)
    for i in range(len(name)):
        p.plot(traj[:,0,i], traj[:,1,i], label=name[i])
    
def printmoon(traj, title_in):
    p.new(title=title_in)
    p.plot(traj[:,0,2] - traj[:,0,1], traj[:,1,2] - traj[:,1,1], label='Moon')

print name

# 3.1 simulate and plot the trajectories of the different particles for one year, dt = 0.0001
traj_euler= simsys(step_eulersym, 0.0001)
printsys(traj_euler, name, 'euler scheme, dt=0.0001a')

# moon in rest frame of earth with different time steps dt
traj_euler2= simsys(step_eulersym, 0.001)
traj_euler3= simsys(step_eulersym, 0.01)

p.new(title='moon in rest frame of earth, different dts')
p.plot(traj_euler[:,0,2] - traj_euler[:,0,1], traj_euler[:,1,2] - traj_euler[:,1,1], label='dt = 0.0001a')
p.plot(traj_euler2[:,0,2] - traj_euler2[:,0,1], traj_euler2[:,1,2] - traj_euler2[:,1,1], label='dt = 0.001a')
p.plot(traj_euler3[:,0,2] - traj_euler3[:,0,1], traj_euler3[:,1,2] - traj_euler3[:,1,1], label='dt = 0.01a')

# 3.2 integrators: simulate and plot trajectories of moon in rest frame of earth for different integrators
traj_euler = simsys(step_euler, 0.01)
traj_eulersym = simsys(step_eulersym, 0.01)
traj_vv = simsysvv(0.01)

p.new(title='moon in rest frame of earth, dt=0.1a')
p.plot(traj_euler[:,0,2] - traj_euler[:,0,1], traj_euler[:,1,2] - traj_euler[:,1,1], label='euler')
p.plot(traj_eulersym[:,0,2] - traj_eulersym[:,0,1], traj_eulersym[:,1,2] - traj_eulersym[:,1,1], label='symplectic euler')
p.plot(traj_vv[:,0,2] - traj_vv[:,0,1], traj_vv[:,1,2] - traj_vv[:,1,1], label='velocity verlet')

# long term stability
traj_euler = simsys(step_euler, 0.01, 10)
traj_eulersym = simsys(step_eulersym, 0.01, 10)
traj_vv = simsysvv(0.01, 10)

p.new(title='moon in rest frame of earth, dt=0.1a, 10 years')
p.plot(traj_euler[:,0,2] - traj_euler[:,0,1], traj_euler[:,1,2] - traj_euler[:,1,1], label='euler')
p.plot(traj_eulersym[:,0,2] - traj_eulersym[:,0,1], traj_eulersym[:,1,2] - traj_eulersym[:,1,1], label='symplectic euler')
p.plot(traj_vv[:,0,2] - traj_vv[:,0,1], traj_vv[:,1,2] - traj_vv[:,1,1], label='velocity verlet')

p.make()

