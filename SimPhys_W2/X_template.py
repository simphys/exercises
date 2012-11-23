#!/usr/bin/python2
# -*- coding:utf-8 -*-

from __future__ import division
import numpy as np

# ---- this section is usefull everytime ----
import os, sys
sys.path.append('../libs')
from evaluation import Plotter
filename = os.path.splitext(os.path.basename(__file__))[0]
p = Plotter(show = True, save = True, pgf = False, name=filename, directory = '')
# -------------------------------------------

# ==== DEFINITIONS ====

# ==== FUNCTIONS ====
    
# ==== CALCULATION ====

# ==== PLOTTING ====
print "test"
p.new(title=u'Title')
p.plot(np.linspace(0,10,100), np.linspace(0,10,100)**2, '+-', label=u'Label % *.2f' % (6,0.42))
p.make()
