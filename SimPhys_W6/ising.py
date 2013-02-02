#!/usr/bin/python2
# -*- coding:utf-8 -*-
# Author: Sebastian Weber, 5 PI Universitaet Stuttgart

from __future__ import division
from libs.simlib import Plotter, Params
import numpy as np
import subprocess
import os.path

# SETUP ###################################################################################################
print '==Setup=='

params = Params()

params.Ising_L = 4
params.Ising_J = 1
params.Ising_H = 0

params.MC_Seed = 42
params.MC_Sweeps = 5000

params.T_Start = 1
params.T_Stop = 5
params.T_StepSize = 0.1

params.Binning_K = 50

useCache = True

p = Plotter(show = True, pdf = False, pgf = False, latex=False, name='ising')

# FUNCTIONS ###############################################################################################

# CALCULATIONS ############################################################################################
print '\n==Calculations=='

fileid = params.id()
filename = './cache/ising_'+params.id()+'.csv'
fileexists = os.path.isfile(filename)

if fileexists and useCache:
    print 'Data file found.'
else:
    print 'Start c++ program ...'
    cmd = './c/Release/SimPhys_W6_C'
    if not os.path.isfile(cmd): cmd = './c/Debug/SimPhys_W6_C'
    if not os.path.isfile(cmd): raise Exception("C++ program not found.")
    args = params.args()
    args +=' -o "'+filename+'"'
    
    proc = subprocess.Popen(cmd+args, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in iter(proc.stdout.readline, ''): print '>> '+line,
    retval = proc.wait()
    print 'Data file created.'

[T,E,M,acceptance] = np.loadtxt(filename, unpack=True)
print 'Data loaded.'

print 'Finished calculations.'

# PLOTS ###################################################################################################
print '\n==Plots=='
p.new(title='Test',xlabel='x',ylabel='y')
p.plot(T,M/params.Ising_L**2)
p.make(ncols=2)
print 'Finished plots.'