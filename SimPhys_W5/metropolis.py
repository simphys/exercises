from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import math

##Functions
def metropolis(N,P,trial_move,phi0,trial_move_param):
    phi=np.zeros(N+1)
    phi[0]=phi0
    actrate=0
    for i in range(N):
        phis=trial_move(phi[i],trial_move_param)
        r=np.random.uniform(0,1)
        if r < min(1,P(phis)/P(phi[i])):
            phi[i+1]=phis
            actrate+=1
        else:
            phi[i+1]=phi[i]
    print actrate/N        
    return phi

def trial_move_x(x,dx):
    r=np.random.uniform(-dx,dx)
    return x+r

def runge(x):
    return 1/(1+x**2)

## Samples for several dx
samples1=metropolis(100000,runge,trial_move_x,0,0.1)
samples2=metropolis(100000,runge,trial_move_x,5,1.0)
samples3=metropolis(100000,runge,trial_move_x,0,10.0)
samples4=metropolis(100000,runge,trial_move_x,0,100.0)
scaling=math.atan(5)-math.atan(-5)
xrunge=np.linspace(-5,5,1000)
yrunge=runge(xrunge)/scaling

plot1=plt.subplot(221)
plot1.hist(samples1,100,range=[-5,5],normed=1)
plot1.plot(xrunge,yrunge,color='red',label='normed runge function')
plt.xlabel('dx=0.1')
plot2=plt.subplot(222)
plot2.hist(samples2,100,range=[-5,5],normed=1)
plot2.plot(xrunge,yrunge,color='red',label='normed runge function')
plt.xlabel('dx=1.0')
plot3=plt.subplot(223)
plot3.hist(samples3,100,range=[-5,5],normed=1)
plot3.plot(xrunge,yrunge,color='red',label='normed runge function')
plt.xlabel('dx=10.0')
plot4=plt.subplot(224)
plot4.hist(samples4,100,range=[-5,5],normed=1)
plot4.plot(xrunge,yrunge,color='red',label='normed runge function')
plt.xlabel('dx=100.0')
plt.show()
