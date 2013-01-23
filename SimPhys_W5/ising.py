#!/usr/bin/python2
# -*- coding:utf-8 -*-
# Author: Sebastian Weber, 5 PI Universitaet Stuttgart

from __future__ import division
from libs.simlib import Plotter
import numpy as np

'''
=== SETUP ===
'''
n = 4
J = 1
H = 0
TStart = 1
TStop = 5
TStep = 0.1
MCSteps = 10000
np.random.seed(42)

'''
=== FUNCTIONS ===
'''
itmp = np.arange(-1,n-1)

#==Exact===
def generateAllStates(n=n):
    nsqrt = n*n
    possibilities = np.mgrid[[slice(-1,3,2) for _ in range(nsqrt)]] # Sorry for the hard to read and memory consuming vectorization, however it makes the calculation fast as hell.
    configurations = possibilities.reshape(n,n,-1).T 
    return configurations

def calcEFromAll(configurations, J=J, H=H, n=n, i=itmp):
    E = -J*np.sum(np.sum(configurations*configurations[:,:,i]+configurations*configurations[:,i,:],axis=1),axis=1)/n**2
    E -= H*np.sum(np.sum(configurations,axis=1),axis=1)/n**2 # Attention: Numpy doesn't like "E /= n**2". It applies integer division!
    return E

def calcMFromAll(configurations, n=n):
    M = np.sum(np.sum(configurations,axis=1),axis=1)/n**2
    return M

def calcMean(O,E,T):
    Exp = np.exp(-E/T)
    PartitionFunction = np.sum(Exp)
    Mean=np.sum(O*Exp)/PartitionFunction
    return Mean

#==MC===
def calcE(configuration, J=J, H=H, n=n, i = itmp):
    E = -J*np.sum(configuration*configuration[:,i]+configuration*configuration[i,:])/n**2
    E -= H*np.sum(configuration)/n**2
    return E

def calcM(configuration, n=n):
    M = np.sum(configuration)/n**2
    return M

def trialFlip(Conf, n=n):
    iRandom = np.random.randint(0,n,2)
    newConf = Conf.copy()
    newConf[iRandom[0],iRandom[1]] *= -1
    return newConf

def metropolisMC(N,trialMove,Phi,T):
    actrate=0
    E = calcE(Phi)
    P = np.exp(-E/T)
    M = calcM(Phi)
    
    arrayE = np.empty(N+1)
    arrayM = np.empty(N+1)
    arrayP = np.empty(N+1)
    arrayE[0] = E
    arrayM[0] = M
    arrayP[0] = P
    
    for i in range(1,N+1):
        newPhi = trialMove(Phi)
        newE = calcE(newPhi) # ToDo: Nur Änderung neu berechnen, es mit zwei Änderungen versuchen
        newP = np.exp(-newE/T)
        r = np.random.uniform(0,1)
        
        if r < min(1,newP/P):
            Phi = newPhi
            E = newE
            M = calcM(Phi)
            P = newP
            actrate+=1
        arrayE[i] = E
        arrayM[i] = M
        arrayP[i] = P

    return arrayE, arrayM, arrayP, actrate/N


'''
=== CALCULATIONS ===
'''
T = np.arange(TStart, TStop+TStep, TStep)


#==Exact===
configurations = generateAllStates(n)
M = calcMFromAll(configurations)
E = calcEFromAll(configurations)

meanE = np.empty_like(T)
meanM = np.empty_like(T)
for i in range(len(T)):
    meanE[i] = calcMean(E,E,T[i])
    meanM[i] = calcMean(M,E,T[i])

#==MC===
configuration = 2*np.random.randint(0,2,(n,n))-np.ones((n,n))

arrayE = np.empty((len(T),MCSteps+1))
arrayM = np.empty((len(T),MCSteps+1))
arrayP = np.empty((len(T),MCSteps+1))
arrayA = np.empty(len(T))

for i in range(len(T)):
    arrayE[i], arrayM[i], arrayP[i], arrayA[i] = metropolisMC(MCSteps,trialFlip,configuration,T[i])
    print T[i]

MC_meanE = np.mean(arrayE,axis=1)
MC_meanM = np.mean(arrayM,axis=1)
MC_acceptance = np.mean(arrayA)

MC_E = arrayE.flat
MC_M = arrayM.flat
MC_P = arrayP

'''
=== PLOTS ===
'''
p = Plotter(show = True, pdf = False, pgf = False, name='ising')

#==Exact===
p.new(name='Mean energy',xlabel='Temperature',ylabel='Energy')
p.plot(T,meanE,label='exact')
p.plot(T,MC_meanE,label='metropolis')
p.new(name='Mean magnetization',xlabel='Temperature',ylabel='Magnetization')
p.plot(T,meanM,label='exact')
p.plot(T,MC_meanM,label='metropolis')
p.new(name='Energy(magnetization)',xlabel='Magnetization',ylabel='Energy')
sort = np.argsort(M)
MC_sort = np.argsort(MC_M)
p.plot(M[sort],E[sort],label='exact')
p.plot(MC_M[MC_sort],MC_E[MC_sort],label='metropolis')
p.new(name='Frequency of probabilities',xlabel='Probability',ylabel='Temperature')
time = (T*np.ones_like(arrayP).T).T
H, xedges, yedges = np.histogram2d(time.flatten(), MC_P.flatten(), bins=(len(T),100))
p.imshow(H, extent=[yedges[0], yedges[-1], xedges[0], xedges[-1]], interpolation='nearest',aspect='auto',origin='lower')

p.make(ncols=2)