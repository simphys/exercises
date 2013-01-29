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
MCSteps = 1000
k = 1
np.random.seed(42)

useExact = True
useMC = True

'''
=== FUNCTIONS ===
'''
itmp = np.arange(-1,n-1)

#==Exact===
def generateAllStates(n=n):
    nsqrt = n*n
    possibilities = np.mgrid[[slice(-1,3,2) for _ in range(nsqrt)]] # Vectorization makes the calculation fast as hell.
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
def metropolisMC(N,Phi,T, n=n, i=itmp):
    # Calculate start values
    E = (-J*np.sum(configuration*configuration[:,i]+configuration*configuration[i,:])-H*np.sum(configuration))/n**2 # Energy
    M = np.sum(Phi)/n**2 # Magnetization
    P = np.exp(-E/T) # Probability
    
    # Initialize variables
    arrayE = np.empty(N+1)
    arrayM = np.empty(N+1)
    arrayP = np.empty(N+1)
    arrayE[0] = E
    arrayM[0] = M
    arrayP[0] = P
    actrate=0
    
    for i in range(1,N+1):
        # Select a spin
        k, l = np.random.randint(0,n,2)
        k_p = k+1
        l_p = l+1
        if k_p >= n: k_p = 0
        if l_p >= n: l_p = 0
        
        # Performe a trial move
        newE = E+2*Phi[k,l]*(J*(Phi[k-1,l]+Phi[k_p,l]+Phi[k,l-1]+Phi[k,l_p])+H)/n**2
        newP = np.exp(-newE/T)
        
        # Attempt to accept the move
        r = np.random.uniform(0,1)
        if r < min(1,newP/P):
            
            # If move accepted, save new values
            E = newE
            P = newP
            M = M-2*Phi[k,l]/n**2
            Phi[k,l] *= -1
            actrate+=1

        # Write current values in array
        arrayE[i] = E
        arrayM[i] = M
        arrayP[i] = P

    return arrayE, arrayM, arrayP, actrate/N

#===ERROR ANALYSIS===
'''
def jackknife(allValues):
    m = allValues.shape[0]
    n = allValues.shape[1]
    meanValues = np.mean(allValues,axis=1)
    
    variances = np.zeros(m)
    onesNxN = np.ones((n,n))
    identNxN = np.eye(n)
    for i in range(m):
        variances[i] = np.mean((np.sum(allValues[i]*onesNxN-meanValues[i]*identNxN,axis=1)/(n-1)-meanValues[i])**2)

    variances *= (n-1)

    return np.sqrt(variances)
'''
'''
def binning(allValues,k = 20):
    [m,n] = allValues.shape[0:2]
    nBlocks=n//k
    
    meanValues = np.mean(allValues[:,:nBlocks*k],axis=1)
    
    allBlocks = allValues[:,:nBlocks*k].reshape((m,k,-1))
    meanBlocks = np.mean(allBlocks, axis = 1)
    
    variances = np.mean((np.subtract(meanBlocks.T,meanValues).T)**2,axis=1)/(nBlocks-1)
    
    return np.sqrt(variances)
    '''

def binning(allValues,k = k):
    nBlocks=len(allValues)//k

    allBlocks = allValues[:nBlocks*k].reshape((k,-1))
    meanBlocks = np.mean(allBlocks,axis=0)
    meanValue = np.mean(meanBlocks)
    
    variance = np.mean((meanBlocks-meanValue)**2)/(nBlocks-1)
    
    return np.sqrt(variance)

def binningAll(arr,k = k):
    n = len(arr)
    result = np.empty(n)
    for i in range(n): result[i] = binning(arr[i],k)
    return result

'''
=== CALCULATIONS ===
'''
T = np.arange(TStart, TStop+TStep, TStep)


#==Exact===
if useExact:
    configurations = generateAllStates(n)
    M = calcMFromAll(configurations)
    E = calcEFromAll(configurations)
    
    meanE = np.empty_like(T)
    meanM = np.empty_like(T)
    meanMabs = np.empty_like(T)
    for i in range(len(T)):
        meanE[i] = calcMean(E,E,T[i])
        meanM[i] = calcMean(M,E,T[i])
        meanMabs[i] = calcMean(abs(M),E,T[i])
    
    print 'Finished exact calculation.'

#==MC===
if useMC:
    configuration = 2*np.random.randint(0,2,(n,n))-np.ones((n,n))
    
    arrayE = np.empty((len(T),MCSteps+1))
    arrayM = np.empty((len(T),MCSteps+1))
    arrayP = np.empty((len(T),MCSteps+1))
    arrayA = np.empty(len(T))
    
    for i in range(len(T)):
        arrayE[i], arrayM[i], arrayP[i], arrayA[i] = metropolisMC(MCSteps,configuration,T[i])
    
    MC_meanE = np.mean(arrayE,axis=1)
    MC_meanM = np.mean(arrayM,axis=1)
    MC_meanMabs = np.mean(abs(arrayM),axis=1)
    MC_acceptance = np.mean(arrayA)
    
    MC_E = arrayE.flat
    MC_M = arrayM.flat
    MC_P = arrayP
    print 'Finished metropolis calculation.'

    MC_errmE = binningAll(arrayE)
    MC_errmM = binningAll(arrayM)
    MC_errmMabs = binningAll(abs(arrayM))
    
    print 'Finished error calculation.'

'''
=== PLOTS ===
'''
p = Plotter(show = True, pdf = False, pgf = False, name='ising')

p.new(name='Mean energy',xlabel='Temperature',ylabel='Energy')
if useExact:
    p.plot(T,meanE,label='exact')
if useMC:
    p.errorbar(T, MC_meanE, yerr=MC_errmE, label='metropolis')

p.new(name='Mean magnetization',xlabel='Temperature',ylabel='Magnetization')
if useExact:
    p.plot(T,meanM,label='exact')
if useMC:
    p.errorbar(T, MC_meanM, yerr=MC_errmM, label='metropolis')

p.new(name='Mean absolute magnetization',xlabel='Temperature',ylabel=r'\vertMagnetization\vert')
if useExact:
    p.plot(T,meanMabs,label='exact')
if useMC:
    p.errorbar(T, MC_meanMabs, yerr=MC_errmMabs, label='metropolis')

p.new(name='Energy(magnetization)',xlabel='Magnetization',ylabel='Energy')
if useExact:
    sort = np.argsort(M)
    p.plot(M[sort],E[sort],label='exact')
if useMC:
    MC_sort = np.argsort(MC_M)
    p.plot(MC_M[MC_sort],MC_E[MC_sort],label='metropolis')

if useMC:
    p.new(name='Frequency of probabilities as a function of Temperature',xlabel='Probability',ylabel='Temperature')
    temp = (T*np.ones_like(arrayP).T).T
    H, xedges, yedges = np.histogram2d(temp.flatten(), arrayP.flatten(), bins=(len(T),100))
    p.imshow(H, extent=[yedges[0], yedges[-1], xedges[0], xedges[-1]], interpolation='nearest',aspect='auto',origin='lower')
    
    """
    p.new(name='Frequency of energies',xlabel='Energy',ylabel='Temperature')
    temp = (T*np.ones_like(arrayE).T).T
    H, xedges, yedges = np.histogram2d(temp.flatten(), arrayE.flatten(), bins=(len(T),100))
    p.imshow(H, extent=[yedges[0], yedges[-1], xedges[0], xedges[-1]], interpolation='nearest',aspect='auto',origin='lower')
    """
    p.new(name='Binning Analysis',xlabel='k',ylabel='error')
    ks = np.arange(1,1000,1)
    error = np.empty((len(ks),len(arrayE)))
    for i in range(len(ks)):
        error[i] = binningAll(arrayE,ks[i])
    p.plot(ks,error)
    

print 'Finished plots.'
p.make(ncols=2)