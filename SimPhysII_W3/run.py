#!/usr/bin/python2
# -*- coding:utf-8 -*-

from __future__ import division
import os, sys, pickle, shutil, fileinput, glob, re, subprocess
import numpy as np
from libeval import Plotter

patternTitle = re.compile(r'@    title "(?P<title>.*)"')
patternComment = re.compile(r'[@#]')

datadict = {}

for item in os.listdir('./'):
    if not item in ['spc','spce','tip3p']:
        continue

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

p = Plotter(show = True, pgf = True, pdf = True, latex=False, sfile=(8,3.5), directory='./report/plots/', loc = 2)

p.new(xlabel=r'time', ylabel=r'number of hydrogen bonds',name='hbnum')
for medium, values in sorted(datadict['hbnum'].iteritems()):
    p.plot(values[1], values[2], '-', label=medium)

p.new(xlabel=r'time', ylabel=r'mean aquare displacement',name='msd')
for medium, values in sorted(datadict['msd'].iteritems()):
    p.plot(values[1], values[2], '-', label=medium)

p.new(xlabel=r'time', ylabel=r'diffusion coefficient',name='diffusion')
for medium, values in sorted(datadict['msd'].iteritems()):
    p.plot(values[1][1:], values[2][1:]/values[1][1:]/6, '-', label=medium)

p.new(xlabel=r'radius', ylabel=r'radial distribution function',name='rdf')
for medium, values in sorted(datadict['rdf'].iteritems()):
    p.plot(values[1], values[2], '-', label=medium)

p.make()