#!/usr/bin/python2
# -*- coding:utf-8 -*-

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

d=np.linspace(5,12,8)
t=[0.148,0.324,0.548,1.052,1.892,3.116,5.116,8.121]
print t[:]

plt.plot(d,t,'+', color='black')
plt.plot(d,t,'-', color='blue')
plt.xlabel('particle number n on one side')
plt.ylabel('timing of the simulation [s]')
plt.show()
