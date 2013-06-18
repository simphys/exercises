#!/usr/bin/python2
# -*- coding:utf-8 -*-

from __future__ import division
import os, sys, pickle, shutil, fileinput, glob, re, subprocess
import numpy as np
from libeval import Plotter, Fitter
from scipy import signal

# === Definitions ===

patternTitle = re.compile(r'@    title "(?P<title>.*)"')
patternComment = re.compile(r'[@#]')
c = ['b','g','r']

# === Parsing ===

datadict = {}

for item in ['spc','spce','tip3p']:

    for task in ['hbnum','msd','rdf']:
        filepath = './'+item+'/'+task+'.xvg'

        if task not in datadict:
            datadict[task] = {}

        datadict[task][item] = [ ]

        with open(filepath, 'r') as f:

            line = f.readline()
            while line != '':

                matchTitle = patternTitle.search(line)
                if matchTitle:
                    title = matchTitle.group("title")

                matchComment = patternComment.search(line)
                if not matchComment:
                    lastPosition = f.tell()-len(line)
                    break

                line = f.readline()

            f.seek(lastPosition)
            table = np.genfromtxt(f, dtype=float, unpack=True, invalid_raise = False)

            datadict[task][item].append(title)
            for col in table:
                datadict[task][item].append(col)

# === Plotting ===

p = Plotter(show = True, pgf = True, pdf = True, latex=False, sfile=(8,3), directory='./report/plots/', loc = 2)
f = Fitter()

# hbnum
print "=== hbnum ==="
p.new(xlabel=r'time [ps]', ylabel=r'number of hydrogen bonds',name='hbnum', loc = 3)
n = 0
for medium, values in sorted(datadict['hbnum'].iteritems()):
    t = values[1]
    num = values[2]

    i = t > 30
    m = [np.mean(num[i]), np.mean(num[i])]
    x = [t[i][0],t[i][-1]]

    print medium+": %f"%m[0]

    p.plot(t, num, c[n]+'-', alpha=0.2, label=medium+', mean: %f'%m[0])
    p.plot(x, m, c[n]+'-',lw=2)
    n+=1

# msd
p.new(xlabel=r'time [ps]', ylabel=r'mean square displacement [nm$^2$]',name='msd', xscale='log',yscale='log')
n = 0
for medium, values in sorted(datadict['msd'].iteritems()):
    p.plot(values[1], values[2], c[n]+'-', label=medium)
    n += 1

# diffusion
print "=== diffusion ==="
p.new(xlabel=r'time [ps]', ylabel=r'diffusion coefficient [10$^{-5}$ cm$^2$/s]',name='diffusion')
n = 0
for medium, values in sorted(datadict['msd'].iteritems()):
    t = values[1][1:]
    D = values[2][1:]/values[1][1:]/6*1e-6*1e4*1e5

    i = np.array(t > 100) & np.array(t < 450)
    m = np.array([np.mean(D[i]), np.mean(D[i])])
    x = [t[i][0],t[i][-1]]

    print medium+": %g"%(m[0])

    p.plot(t, D, c[n]+'-', label=medium+', mean: %g'%(m[0]))
    p.plot(x, m, c[n]+'-',lw=2)
    n += 1

# rdf
print "=== rdf ==="
p.new(xlabel=r'distance [nm]', ylabel=r'radial distribution function',name='rdf')
n = 0
for medium, values in sorted(datadict['rdf'].iteritems()):
    i1 = signal.find_peaks_cwt(values[2], np.arange(1,20),noise_perc=1)[0]
    i2 = signal.find_peaks_cwt(values[2], np.arange(20,40),noise_perc=1)[1:3]
    peakindex = np.append(i1,i2)
    mx = values[1][peakindex]
    my = values[2][peakindex]
    print medium+":"
    for i in range(0,3):
        print "(%f | %f)"%(mx[i],my[i])
    p.plot(mx, my, c[n]+'*',ms=11, mew=1)
    p.plot(values[1], values[2], c[n]+'-', label=medium)
    n += 1

p.make()




'''f.loadFunction(lambda x, *p: p[1]*x**p[0], [1,1])
p.new(xlabel=r'time [ps]', ylabel=r'mean aquare displacement',name='regimes', xscale='log',yscale='log')
for medium, values in sorted(datadict['msd'].iteritems()):
    i = values[1] < 600

    j = values[1][i] < 0.3
    f.loadData(values[1][i][j], values[2][i][j], scale='log')
    x = np.linspace(0.1, 10,2)
    p.plot(x, f(x), 'k-')

    j = values[1][i] > 0.5
    f.loadData(values[1][i][j], values[2][i][j], scale='log')
    x = np.linspace(0.1, 1,2)
    p.plot(x, f(x), 'k-')

    p.plot(values[1][i], values[2][i], '-', label=medium)'''