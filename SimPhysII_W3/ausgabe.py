#!/usr/bin/python2
# -*- coding:utf-8 -*-

from __future__ import division
import os, sys, pickle, shutil, fileinput, glob, re, subprocess
import numpy as np
from libeval import Plotter, Fitter

patternTitle = re.compile(r'@    title "(?P<title>.*)"')
patternComment = re.compile(r'[@#]')

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

p = Plotter(show = True, pgf = False, pdf = False, latex=False, sfile=(8,3.5), directory='./report/plots/', loc = 2)
f = Fitter()

p.new(xlabel=r'time [ps]', ylabel=r'number of hydrogen bonds',name='hbnum')
for medium, values in sorted(datadict['hbnum'].iteritems()):
    t = values[1]
    num = values[2]

    i = t > 100
    m = [np.mean(num[i]), np.mean(num[i])]
    x = [t[i][0],t[i][-1]]

    p.plot(t, num, '-', label=medium+', mean: %f'%m[0])
    p.plot(x, m, 'k-')

p.new(xlabel=r'time [ps]', ylabel=r'mean aquare displacement',name='msd')
for medium, values in sorted(datadict['msd'].iteritems()):
    p.plot(values[1], values[2], '-', label=medium)

f.loadFunction(lambda x, *p: p[1]*x**p[0], [1,1])
p.new(xlabel=r'time [ps]', ylabel=r'mean aquare displacement',name='regimes', xscale='log',yscale='log')
for medium, values in sorted(datadict['msd'].iteritems()):
    i = values[1] < 600

    '''j = values[1][i] < 0.3
    f.loadData(values[1][i][j], values[2][i][j], scale='log')
    x = np.linspace(0.1, 10,2)
    p.plot(x, f(x), 'k-')

    j = values[1][i] > 0.5
    f.loadData(values[1][i][j], values[2][i][j], scale='log')
    x = np.linspace(0.1, 1,2)
    p.plot(x, f(x), 'k-')'''

    p.plot(values[1][i], values[2][i], '-', label=medium)

p.new(xlabel=r'time [ps]', ylabel=r'diffusion coefficient',name='diffusion')
for medium, values in sorted(datadict['msd'].iteritems()):
    t = values[1][1:]
    D = values[2][1:]/values[1][1:]/6

    i = t > 200
    m = [np.mean(D[i]), np.mean(D[i])]
    x = [t[i][0],t[i][-1]]

    p.plot(t, D, '-', label=medium+', mean: %f'%m[0])
    p.plot(x, m, 'k-')

p.new(xlabel=r'radius', ylabel=r'radial distribution function',name='rdf')
for medium, values in sorted(datadict['rdf'].iteritems()):
    p.plot(values[1], values[2], '-', label=medium)

p.make()