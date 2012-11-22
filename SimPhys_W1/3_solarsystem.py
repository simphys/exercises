#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import matplotlib.pyplot as p
import os

# ~~~~~~~ this section is usefull whatever you do... ~~~~~~~~~~~
plotid = 0
def plotStart(xlabel=u'', ylabel=u'', xscale = 'linear', yscale = 'linear'):
    global plotid
    p.figure(plotid,figsize=(10,4))
    p.xlabel(xlabel)
    p.ylabel(ylabel)
    p.xscale(xscale)
    p.yscale(yscale)
    p.grid()
    plotid += 1
def plotEnd(show = True, save = True, pgf = False, directory = '', name = ''):
    global plotid
    p.rc('legend', fontsize='small')
    p.legend(shadow=0, loc='best')
    if save:
        if not directory: directory = './plots/'
        if not name: name = os.path.splitext(os.path.basename(__file__))[0]+'_%0*i'%(2,plotid)
        if not os.path.isdir(directory): os.mkdir(directory)
        if pgf: p.savefig(directory+name+'.pgf')
        else: p.savefig(directory+name+'.pdf', bbox_inches='tight')
    if not show: p.close()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

name, x_init, v_init, m, g = np.load('solar_system.npy')

#print name
#print x_init
#print v_init
#print m
#print g

def compute_forces(x, m):
    global g
    l, c = x.shape
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
    plotStart()
    for i in range(len(name)):
        p.plot(traj[:,0,i], traj[:,1,i], label=name[i])
    p.title(title_in)
    plotEnd()
    
def printmoon(traj, title_in):
    plotStart()
    p.plot(traj[:,0,2] - traj[:,0,1], traj[:,1,2] - traj[:,1,1], label='Moon')
    p.title(title_in)
    plotEnd()

print name

# 3.1 simulate and plot the trajectories of the different particles for one year, dt = 0.0001
traj_euler= simsys(step_eulersym, 0.0001)
printsys(traj_euler, name, 'euler scheme, dt=0.0001a')

# moon in rest frame of earth with different time steps dt
traj_euler2= simsys(step_eulersym, 0.001)
traj_euler3= simsys(step_eulersym, 0.01)

plotStart()
p.plot(traj_euler[:,0,2] - traj_euler[:,0,1], traj_euler[:,1,2] - traj_euler[:,1,1], label='dt = 0.0001a')
p.plot(traj_euler2[:,0,2] - traj_euler2[:,0,1], traj_euler2[:,1,2] - traj_euler2[:,1,1], label='dt = 0.001a')
p.plot(traj_euler3[:,0,2] - traj_euler3[:,0,1], traj_euler3[:,1,2] - traj_euler3[:,1,1], label='dt = 0.01a')
p.title('moon in rest frame of earth, different dts')
plotEnd()

# 3.2 integrators: simulate and plot trajectories of moon in rest frame of earth for different integrators
traj_euler = simsys(step_euler, 0.01)
traj_eulersym = simsys(step_eulersym, 0.01)
traj_vv = simsysvv(0.01)

plotStart()
p.plot(traj_euler[:,0,2] - traj_euler[:,0,1], traj_euler[:,1,2] - traj_euler[:,1,1], label='euler')
p.plot(traj_eulersym[:,0,2] - traj_eulersym[:,0,1], traj_eulersym[:,1,2] - traj_eulersym[:,1,1], label='symplectic euler')
p.plot(traj_vv[:,0,2] - traj_vv[:,0,1], traj_vv[:,1,2] - traj_vv[:,1,1], label='velocity verlet')
p.title('moon in rest frame of earth, dt=0.1a')
plotEnd()

# long term stability
traj_euler = simsys(step_euler, 0.01, 10)
traj_eulersym = simsys(step_eulersym, 0.01, 10)
traj_vv = simsysvv(0.01, 10)

plotStart()
p.plot(traj_euler[:,0,2] - traj_euler[:,0,1], traj_euler[:,1,2] - traj_euler[:,1,1], label='euler')
p.plot(traj_eulersym[:,0,2] - traj_eulersym[:,0,1], traj_eulersym[:,1,2] - traj_eulersym[:,1,1], label='symplectic euler')
p.plot(traj_vv[:,0,2] - traj_vv[:,0,1], traj_vv[:,1,2] - traj_vv[:,1,1], label='velocity verlet')
p.title('moon in rest frame of earth, dt=0.1a, 10 years')
plotEnd()

p.show()

