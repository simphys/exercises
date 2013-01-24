from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import math

## Functions
## Calculating f(x)
def runge(x):
    return 1/(1+x**2)

## Exact integral of f(x) 
def runge_int_ex(a,b):
    return math.atan(b)-math.atan(a)

## Simple sampling
def simple_sampling(f,a,b,N):
    xi=np.random.uniform(a,b,N)
    return (b-a)*np.sum(f(xi))/N, (b-a)*np.std(f(xi))/math.sqrt(N)

## Plot of the runge function
#x=np.linspace(-5,5,1000)
#y=runge(x)
#plt.plot(x,y)
#plt.xlabel('x')
#plt.ylabel('f(x)')
#plt.xlim([-5,5])
#plt.ylim([0,1.1])
#plt.show()


exact_solution=runge_int_ex(-5,5)
print exact_solution

## Task on simple sampling

#print simple_sampling(runge,-5,5,1000)

erg_I=[]
for i in range(2,21):
    erg_I.append(simple_sampling(runge,-5,5,2**i))

erg_I=np.array(erg_I)
act_error=np.abs(erg_I[:,0]-exact_solution)
N=2**np.arange(2,21)
plt.semilogx(N,act_error,label="actual error")
plt.semilogx(N,erg_I[:,1],label="statistical error")
plt.xlabel('Number of integrationsteps')
plt.ylabel('Error')
plt.legend()
plt.show()
