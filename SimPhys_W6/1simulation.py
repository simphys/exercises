#!/usr/bin/python2
# -*- coding:utf-8 -*-
# Author: Sebastian Weber, 5 PI Universitaet Stuttgart

from __future__ import division
from libs.simlib import Plotter, Params
import numpy as np
import subprocess
import os.path
from matplotlib.colors import LogNorm

# SETUP ###################################################################################################
print '==Setup=='

params = Params()

params.Ising_J = 1
params.Ising_H = 0

params.MC_Seed = 42
params.MC_Sweeps = 10000

params.T_Start = 1
params.T_Stop = 5
params.T_StepSize = 0.1

Ising_L = [4,16,64]
Binning_K = [50,200,800]

useCache = True

p = Plotter(show = True, pdf = False, pgf = False, latex=False, name='1simulation')

# FUNCTIONS ###############################################################################################
#==Exact===
itmp = np.arange(-1,4-1)
def generateAllStates(n=4):
    nsqrt = n*n
    possibilities = np.mgrid[[slice(-1,3,2) for _ in range(nsqrt)]] # Vectorization makes the calculation fast as hell.
    configurations = possibilities.reshape(n,n,-1).T 
    return configurations

def calcEFromAll(configurations, J=params.Ising_J, H=params.Ising_H, n=4, i=itmp):
    E = -J*np.sum(np.sum(configurations*configurations[:,:,i]+configurations*configurations[:,i,:],axis=1),axis=1)
    E -= H*np.sum(np.sum(configurations,axis=1),axis=1) # Attention: Numpy doesn't like "E /= n**2". It applies integer division!
    return E

def calcMFromAll(configurations, n=4):
    M = np.sum(np.sum(configurations,axis=1),axis=1)
    return M

def calcMean(O,E,T,n=4):
    Exp = np.exp(-E/T)
    PartitionFunction = np.sum(Exp)
    Exp = np.exp(-E/T)
    Mean=np.sum(O*Exp)/PartitionFunction/n**2
    return Mean

#===Error analysis===
def binning(allValues,k):
    nBlocks=len(allValues)//k

    allBlocks = allValues[:nBlocks*k].reshape((-1,k))
    meanBlocks = np.mean(allBlocks,axis=1)
    meanValue = np.mean(meanBlocks)
    
    variance = np.mean((meanBlocks-meanValue)**2)/(nBlocks-1)
    
    return np.sqrt(variance)

def errorAll(func,arr,k):
    n = len(arr)
    result = np.empty(n)
    for i in range(n): result[i] = func(arr[i],k)
    return result

#===Others===
def load(Ising_L):
    params.Ising_L = Ising_L
    filename = './cache/ising_'+params.id()+'.csv'
    fileexists = os.path.isfile(filename)
    
    if fileexists and useCache:
        print ' Data file found.'
    else:
        print ' Run c++ program ...'
        cmd = './c/Release/SimPhys_W6_C'
        if not os.path.isfile(cmd): cmd = './c/Debug/SimPhys_W6_C'
        if not os.path.isfile(cmd): raise Exception("C++ program not found.")
        args = params.args()
        args +=' -o "'+filename+'"'
        
        proc = subprocess.Popen(cmd+args, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in iter(proc.stdout.readline, ''): print ' >> '+line,
        proc.wait()
        print ' Data file created.'
    print ' Data loaded.'
    
    return np.loadtxt(filename);

def unique(arr):
    order = np.lexsort(arr.T)
    arr = arr[order]
    diff = np.diff(arr, axis=0)
    ui = np.ones(len(arr), 'bool')
    ui[1:] = (diff != 0).any(axis=1)
    return arr[ui]

# CALCULATIONS ############################################################################################
print '\n==Calculations=='

#===Exact===
T_exact = np.arange(params.T_Start, params.T_Stop+params.T_StepSize, params.T_StepSize)
configurations = generateAllStates(4)

allM_exact = calcMFromAll(configurations)
allE_exact = calcEFromAll(configurations)

E_exact = np.empty_like(T_exact)
M_exact = np.empty_like(T_exact)
for i in range(len(T_exact)):
    E_exact[i] = calcMean(allE_exact,allE_exact,T_exact[i])
    M_exact[i] = calcMean(abs(allM_exact),allE_exact,T_exact[i])

print 'Finished exact calculation.'
    
#===Metropolis===
T = []
A = []
arrE = []
arrM = []
E = []
M = []
errE = []
errM = []

for i in range(len(Ising_L)):
    arr = load(Ising_L[i])
    V = Ising_L[i]**2 
    T.append(arr[::2,0])
    A.append(arr[::2,1])
    arrE.append(arr[::2,2:]/V)
    arrM.append(arr[1::2,2:]/V)
    E.append(np.mean(arrE[i],axis=1))
    M.append(np.mean(arrM[i],axis=1))
    errE.append(errorAll(binning,arrE[i],Binning_K[i]))
    errM.append(errorAll(binning,arrM[i],Binning_K[i]))

print 'Finished calculations.'

# PLOTS ###################################################################################################
print '\n==Plots=='

p.new(title='Mean energy',xlabel='Temperature',ylabel='Energy')
p.plot(T_exact,E_exact,label='L = 4, exact')
for i in range(len(Ising_L)): p.errorbar(T[i], E[i], yerr=errE[i], label='L=%s, MC'%Ising_L[i])

p.new(title='Mean absolute magnetization',xlabel='Temperature',ylabel=r'$\vert Magnetization \vert$')
p.plot(T_exact,M_exact,label='L = 4, exact')
for i in range(len(Ising_L)): p.errorbar(T[i], M[i], yerr=errM[i], label='L=%s, MC'%Ising_L[i])
#p.make(ncols=1,show=False)

for i in range(len(Ising_L)):
    p.new(title='Frequency of energies (L=%s)'%Ising_L[i],xlabel='Energy',ylabel='Temperature')
    temp = (T[i]*np.ones_like(arrE[i]).T).T
    H, xedges, yedges = np.histogram2d(temp.flatten(), arrE[i].flatten(), bins=(len(T[i]),100))
    p.imshow(H, extent=[yedges[0], yedges[-1], xedges[0], xedges[-1]], interpolation='nearest',aspect='auto',origin='lower',norm=LogNorm())
    
    p.new(title='Frequency of magnetizations (L=%s)'%Ising_L[i],xlabel='Absolute magnetization',ylabel='Temperature')
    temp = (T[i]*np.ones_like(arrM[i]).T).T
    H, xedges, yedges = np.histogram2d(temp.flatten(), arrM[i].flatten(), bins=(len(T[i]),100))
    p.imshow(H, extent=[yedges[0], yedges[-1], xedges[0], xedges[-1]], interpolation='nearest',aspect='auto',origin='lower',norm=LogNorm())
#p.make(ncols=2,show=False)

for i in range(len(Ising_L)):
    p.new(title='Binning error (L=%s)'%Ising_L[i],xlabel='k',ylabel='error')
    ks = np.arange(1,800,1)
    error = np.empty((len(ks),len(arrE[i])))
    for j in range(len(ks)):
        error[j] = errorAll(binning,arrE[i],ks[j])
    p.plot(ks,error)
    
    p.new(title='Energy of absolute magnetization (L=%s)'%Ising_L[i],xlabel='Absolute magnetization',ylabel='Energy')
    [X, Y] = unique(np.array([np.round(arrM[i].ravel(),2),np.round(arrE[i].ravel(),2)]).T).T
    p.plot(X,Y,'o',alpha=0.1)
p.make(ncols=2,show=True)

print 'Finished plots.'