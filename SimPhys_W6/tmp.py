#!/usr/bin/python2
# -*- coding: utf-8 -*-

# ---- Autor ----
# Sebastian Weber

# Schwingungsgleichungen
# --- mechanisch ---
# m*x`` + b*x` + D*x = F(t)
# <=> x`` + 2*g*x` + w0^2*x = F(t)/m
# z.B. mit F(t)/m = alpha*sin(w*t)
# --- elektromagnetisch ---
# L*q`` + R*q` + 1/C*q = U(t)

from __future__ import division
import numpy as np
from libeval import Plotter, Fitter

# === Constants ===
# Charakterstika der Schwingung
m=1
b=0.1
D=1
Fmax = 1

F = lambda t: Fmax*np.sin(1*t)

# Konstanten für die Simulation
tend = 100
dt = 0.1
accuracy = 0.001

# Anfangswerte
t = 0
x = 0
v = 1
a = (-D*x)/m

# === Functions ===
def step_vv(x, v, a, t, dt=dt):
    global m, b, D, F
    # update positions
    x += v*dt + 0.5*a * dt*dt
    # half update of the velocity
    v += 0.5*a * dt
    # compute new acceleration
    a = (-D*x-b*v+F(t))/m
    # second half update of the velocity
    v += 0.5*a * dt
    return x, v, a

def simulate(x,v,a,t):
    steps = tend//dt
    array_t = np.empty(steps)
    array_x = np.empty(steps)
    array_v = np.empty(steps)
    array_a = np.empty(steps)
    array_E = np.empty(steps)
    
    for n in range(len(array_t)):
        x,v,a = step_vv(x, v, a,t)
        t += dt
        array_t[n]=t
        array_x[n]=x
        array_v[n]=v
        array_a[n]=a
        array_E[n]=0.5*m*v*v+0.5*D*x*x
        
    return array_t, array_x, array_v, array_a, array_E

def getAmplitude(x,v,a,t, acc = accuracy):
    x_max = 0
    x_min = 0
    x_pre = x
    amp = 0
    amp_pre = 0
    nulls = 0
    
    for _n in range(100000):
        x,v,a = step_vv(x, v, a,t)
        t += dt
        x_max = max(x_max,x)
        x_min = min(x_min,x)
        
        # Nulldurchgang registieren
        if x_pre*x <= 0: nulls += 1
        x_pre = x
        
        # Nach einer Periode (3 Nulldurchgänge), Amplitude auswerten
        if nulls > 2:
            amp = (x_max - x_min)/2
            if abs(amp-amp_pre) < acc: break
            amp_pre = amp
            x_max = 0
            x_min = 0
            nulls = 1
            
    return amp
    
# === Calculation ===
x_s = []
v_s = []
a_s = []
E_s = []

for w in np.arange(0,2,0.5):
    F = lambda t: Fmax*np.sin(w*t)
    array_t, array_x, array_v, array_a, array_E = simulate(x,v,a,t)
    amp = getAmplitude(x,v,a,t)
    x_s.append([array_t, array_x, r'$\omega = %.2f$'%w])
    v_s.append([array_t, array_v, r'$\omega = %.2f$'%w])
    a_s.append([array_t, array_a, r'$\omega = %.2f$'%w])
    E_s.append([array_t, array_E, r'$\omega = %.2f$'%w])

w_s = np.arange(0.001,2,0.1)
amp_s = []

for w in w_s:
    F = lambda t: Fmax*np.sin(w*t)
    amp = getAmplitude(x,v,a,t)
    amp_s.append(amp)

# === Plotting ===
p = Plotter(show = True, pgf = False, pdf = False, loc=1)

# --- Auslenkung ---
p.new(title='Auslenkung',xlabel='t',ylabel='x')
for v,w,name in x_s: p.plot(v,w,label=name)

# --- Auslenkung ---
p.new(title='Plots für ausgeschaltete Anregung',xlabel='t')
v,w,name = x_s[0]
p.plot(v,w,label='Auslenkung x')
v,w,name = v_s[0]
p.plot(v,w,label='Geschwindigkeit v')
v,w,name = a_s[0]
p.plot(v,w,label='Beschleunigung a')

# --- Energie ---
p.new(title='Gesamte Energie', yscale='log',xlabel='t',ylabel='E')
for v,w,name in E_s: p.plot(v,w,label=name)

# --- Amplitude ---
p.new(title='Amplitude',xlabel=r'$\omega$',ylabel=r'$\hat{x}$')
f_amp = lambda x, *p: p[0]/np.sqrt((p[1]*p[1]-x*x)**2+(2*p[2]*x)**2)
alpha = Fmax/m
omeganull = np.sqrt(D/m)
gamma = b/m/2
f = Fitter(f_amp, [alpha, omeganull, gamma])
f.loadData(w_s, amp_s, scale='linear')
p.plot(w_s,amp_s,'+',label=r'Simulation mit $\alpha = %.4f$, $\gamma = %.4f$, $\omega_0 = %.4f$'%(alpha, gamma, omeganull))
p.plot(f.x,f.y,'-',label=r'Fitt mit $\alpha = %.4f \pm %.4f$, $\gamma = %.4f \pm %.4f$, $\omega_0 = %.4f \pm %.4f$'%(f.params[0],f.std[0],f.params[2],f.std[2],f.params[1],f.std[1]))

p.make(ncols = 2, savewindow = False)