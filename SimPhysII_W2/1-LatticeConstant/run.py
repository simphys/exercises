#!/usr/bin/python2
# -*- coding:utf-8 -*-

from __future__ import division
import os, sys, pickle, shutil, fileinput, glob, re, subprocess
import numpy as np
from libeval import Plotter

# ==== DEFINE ====
GenerateNewConf = True
GenerateNewData = True
RunCalculations = True

latticeConstants =  {'LDA': [4.95,5.05,5.15,5.25,5.0,5.1,5.2,5.3,5.4,5.43,5.5,5.6,5.7,5.8,5.9],
    'GGA': [5.35,5.25,5.0,5.1,5.2,5.3,5.4,5.43,5.5,5.6,5.7,5.8,5.9]}


datadict = {'LDA': [np.array([]),np.array([]),np.array([])],
    'GGA' : [np.array([]),np.array([]),np.array([])]}
patternE = re.compile(r"(?P<txt>Total =.*?)(?P<energy>[-+]?[0-9]*\.?[0-9]+)")
patternV = re.compile(r"(?P<txt>Cell volume =.*?)(?P<volume>[-+]?[0-9]*\.?[0-9]+)")
patternM = re.compile(r"(?P<txt>Begin CG move =.*?)(?P<move>[0-9]+)")
patternL = re.compile(r"(?P<txt>LatticeConstant.*?)(?P<lattice>[-+]?[0-9]*\.?[0-9]+)")

p = Plotter(show = True, pgf = True, pdf = False, latex=False, sfile=(8,3.5), directory='../report/plots/', loc = 4)

# ==== SETUP ====
if GenerateNewConf:
    for item in os.listdir('./'):
        if not item in ['LDA','GGA']:
            continue

        for l in latticeConstants[item]:
            path = './'+item+'/%.3s'%(l*100)
            if not os.path.isdir(path):
                shutil.copytree('./'+item+'/raw',path)
                fdf = glob.glob(path+'/*.fdf')[0]
                for line in fileinput.FileInput(fdf,inplace=1):
                    line = re.sub(patternL,r'\g<txt>%.2f'%l,line,count=1)
                    print line,

# ==== CALC ====
runnumber = 0
if GenerateNewData:
    for root, _dirs, filenames in os.walk('./'):
        v = ''
        e = ''
        l = ''

        try:
            type, number = root.split('/')[1:]
            fdf = glob.glob(root+"/*.fdf")[0]
        except:
            continue

        if not (type in ['LDA','GGA']) or not number.isdigit():
            continue

        runnumber += 1

        try:
            calculated = glob.glob(root+"/*.xyz")
        except:
            continue

        if not calculated and RunCalculations:
            head, tail = os.path.split(fdf)
            cmd = 'mpirun -np 2 siesta < '+tail+' | tee OUT.out'
            proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,cwd=root)
            for line in iter(proc.stdout.readline, ''):
                matchM = patternM.search(line)
                if matchM:
                    print "%s, run %2s) Begin CG move = %2s"%(type, runnumber,matchM.group("move"))
            proc.wait()

        try:
            out = glob.glob(root+"/*.out")[0]
        except:
            continue

        with open(out, 'r') as f:
            file = f.read()
            matchV = patternV.search(file)
            matchE = patternE.search(file)
            matchL = patternL.search(file)
            if matchV and matchE and matchL:
                v = matchV.group("volume")
                e = matchE.group("energy")
                l = matchL.group("lattice")

        try:
            volume = float(v)
            energy = float(e)
            lattice = float(l)
        except:
            continue

        datadict[type][0] = np.append(datadict[type][0], volume)
        datadict[type][1] = np.append(datadict[type][1], energy)
        datadict[type][2] = np.append(datadict[type][2], lattice)

        print "%s, run %2s) lattice constant = %f, cell volume = %f, total energy = %f \n"%(type,
            runnumber,lattice,volume,energy)

    with open('data.pkl', 'wb') as f:
        pickle.dump(datadict, f)

else:
    with open('data.pkl', 'rb') as f:
        datadict = pickle.load(f)

# ==== PLOT ====
for type, [volume,energy,lattice] in datadict.items():
    sort = np.argsort(volume)
    lattice = np.array(lattice)[sort]
    volume = np.array(volume)[sort]
    energy = np.array(energy)[sort]
    p.new(xlabel=r'cell volume [Ang$^3$]', ylabel=r'total energy [eV]',name=type+"_energy")
    p.plot(volume, energy, '-', label=type)
    p.new(xlabel=r'cell volume [Ang$^3$]',ylabel=r'lattice constant [Ang]',name=type+"_lattice")
    p.plot(volume, lattice, 'x', label=type)

p.make(ncols = 2)