#!/usr/bin/python2
# -*- coding:utf-8 -*-

from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('GtkAgg')
import matplotlib.pyplot as p
import os

# ~~~~~~~ this section is usefull whatever you do... ~~~~~~~~~~~
plotid = 0
def plotStart(xlabel=u'x in m', ylabel=u'y in m', xscale = 'linear', yscale = 'linear'):
    global plotid
    p.figure(plotid,figsize=(10,4))
    p.xlabel(xlabel)
    p.ylabel(ylabel)
    p.xscale(xscale)
    p.yscale(yscale)
    p.grid()
    plotid += 1
def plotEnd(show = False, save = True, pgf = True, directory = '', name = ''):
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

# ==== DEFINITIONS ====
m = 2.
g = 9.81
dt = 0.1

# ==== FUNCTIONS ====
def compute_forces(x,v,vw = np.array([0.,0.]), gamma = 0.):
    global m, g
    return np.array([0.,-m*g])-gamma*(v-vw)

def step_euler(x,v,f):
    global m, dt
    v = v + f/m *dt
    x = x + v*dt
    return x, v

def calc(vw = np.array([0.,0.]), gamma = 0.):
    global dt
    traj = np.array([[0.,0.]])
    vel = np.array([[50.,50.]])
    t = 0
    while (traj[-1,1] >= 0):
        x, v = step_euler(traj[-1],vel[-1],compute_forces(traj[-1],vel[-1],vw, gamma))
        traj = np.vstack((traj,x))
        vel = np.vstack((traj,v))
        t += dt
    return traj, vel, t
    
# ==== CALCULATION ====
traj_no_fric, vel_no_fric, t_no_fric = calc([0.,0.], 0.)
traj_with_fric_vw0, vel_with_fric_vw0, t_with_fric_vw0 = calc([0.,0.], 0.1)
traj_with_fric_vw50, vel_with_fric_vw50, t_with_fric_vw50 = calc([-50.,0.], 0.1)

vwvarious = []
for vw in np.linspace(-300,10,10): vwvarious.append([vw, calc([vw,0.], 0.1)[0]])

# ==== PLOTTING ====
plotStart()
p.title(u'Trajectories of a cannonball')
p.plot(traj_no_fric[:,0], traj_no_fric[:,1], 'r+-', label=u'no friction \n flying time: %f s' % (t_no_fric))
p.plot(traj_with_fric_vw0[:,0], traj_with_fric_vw0[:,1], 'b+-', label=u'with friction, no wind \n flying time: %f s' % (t_with_fric_vw0))
p.plot(traj_with_fric_vw50[:,0], traj_with_fric_vw50[:,1], 'g+-', label=u'with friction, strong wind \n flying time: %f s' % (t_with_fric_vw50))
plotEnd()

plotStart()
p.title(u'Trajectories of a cannonball at various wind speeds vw')
for vw,traj in vwvarious: p.plot(traj[:,0], traj[:,1], '+-', label=u'vw =% *.2f' % (6,vw))
plotEnd()

p.show()