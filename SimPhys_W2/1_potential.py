#!/usr/bin/python2
# -*- coding:utf-8 -*-

from __future__ import division
import numpy as np
from libs.simlib import Plotter

p = Plotter(show = True, save = False, pgf = True, name='1_potential', directory = '')

# ==== DEFINITIONS ====
EPS = 1
SIG = 1
NParticles = 1000

# ==== FUNCTIONS ====
def compute_lj_potential(rij, eps = EPS, sig = SIG):
    q = sig/np.linalg.norm(rij)
    return 4*eps*(q**12-q**6)

def compute_lj_force(rij, eps = EPS, sig = SIG):
    norm = np.linalg.norm(rij)
    q = sig/norm
    return 4*eps*(12*q**11-6*q**5)*q/norm**2*rij

# ==== CALCULATION ====
d = np.zeros((NParticles,3))
d[:,0] = np.linspace(0.85,2.5,NParticles)
potential = np.array(map(compute_lj_potential, d))
force = np.array(map(compute_lj_force, d))

# ==== PLOTTING ====
p.new(title=u'LJ potential', xlabel ='distance', ylabel = 'potential')
p.plot(d[:,0], potential, '-', label=u'potential')
p.new(title=u'LJ force', xlabel ='distance', ylabel = '1st component of force')
p.plot(d[:,0], force[:,0], '-', label=u'force')
p.make()