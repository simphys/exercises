#!/usr/bin/python2
# -*- coding:utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import matplotlib
matplotlib.use('GtkAgg')
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['text.latex.unicode']=True
import matplotlib.pyplot as p
import numpy as np
import os

class Plot(object):
    def __init__(self,title,name,xlabel,ylabel,xscale,yscale,aspect,show,pdf,pgf,directory):
        self.curves = []
        self.title = title
        self.name = name
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.xscale = xscale
        self.yscale = yscale
        self.aspect = aspect
        self.show = show
        self.pdf = pdf
        self.pgf = pgf
        self.dir = directory
        self.legend = False
        
    def addCurve(self,*args, **kwargs):
        self.curves.append([args,kwargs])
        if kwargs.get('label',None): self.legend = True

class Plotter(object):
    def __init__(self, **kwargs):
        self.__global_xlabel  = kwargs.get('xlabel','')
        self.__global_ylabel  = kwargs.get('ylabel','')
        self.__global_xscale  = kwargs.get('xscale','linear')
        self.__global_yscale  = kwargs.get('yscale','linear')
        self.__global_aspect  = kwargs.get('aspect','auto')
        self.__global_show    = kwargs.get('show',True)
        self.__global_pdf     = kwargs.get('pdf',False)
        self.__global_pgf     = kwargs.get('pgf',False)
        self.__global_direc   = kwargs.get('directory','./plots/')
        self.__global_name    = kwargs.get('name','plot')
        
        self.__reset()
    
    def new(self, **kwargs):
        xlabel  = kwargs.get('xlabel',self.__global_xlabel)
        ylabel  = kwargs.get('ylabel',self.__global_ylabel)
        xscale  = kwargs.get('xscale',self.__global_xscale)
        yscale  = kwargs.get('yscale',self.__global_yscale)
        aspect  = kwargs.get('aspect',self.__global_aspect)
        show    = kwargs.get('show',self.__global_show)
        pdf     = kwargs.get('pdf',self.__global_pdf)
        pgf     = kwargs.get('pgf',self.__global_pgf)
        direc   = kwargs.get('directory',self.__global_direc)
        name    = kwargs.get('name',self.__global_name+'_%0*i'%(2,self.__nr_id))
        title   = kwargs.get('title',name)
        
        self.__plots.append(Plot(r'\verb#Plot %i) '%(self.__nr_id)+title+r'#',name,xlabel,ylabel,xscale,yscale,aspect,show,pdf,pgf,direc))
        
        if show: self.__nr_show += 1
        self.__nr_id += 1
    
    def plot(self, *args, **kwargs):
        self.__plots[-1].addCurve(*args, **kwargs)
    
    def make(self, ncols = 2, swindow = (15,10), sfile = (8,3.5)):        
        for n,plot in enumerate(self.__plots):
            if plot.pdf or plot.pgf: self.__save(n,plot,sfile)
        
        if self.__nr_show > 0:
            if self.__nr_show == 1: f, axarr = p.subplots(1, 1, figsize=swindow)
            else: f, axarr = p.subplots(int(np.ceil(1.*self.__nr_show/ncols)), ncols, figsize=swindow)
            nplots = 0
            for n,plot in enumerate(self.__plots):
                if plot.show:
                    self.__show(n,plot,np.ravel(axarr)[nplots])
                    nplots += 1
            f.tight_layout()
            p.show()
        
        self.__reset()
    
    def __show(self,n,plot, ax):
        ax.set_title(plot.title,fontsize='medium', fontweight='bold', x=.05, y =1., ha = 'left', va='bottom')
        ax.set_xlabel(plot.xlabel,fontsize='small')
        ax.set_ylabel(plot.ylabel,fontsize='small')
        ax.tick_params(labelsize='small')
        ax.set_xscale(plot.xscale)
        ax.set_yscale(plot.yscale)
        ax.grid()
        for args, kwargs in plot.curves: ax.plot(*args, **kwargs)
        ax.set_aspect(plot.aspect)
        if plot.legend: ax.legend(shadow=0, loc='best',fontsize='small')
        
    def __save(self,n,plot,sfile):
        p.figure(figsize=sfile)
        
        p.xlabel(plot.xlabel)
        p.ylabel(plot.ylabel)
        p.xscale(plot.xscale)
        p.yscale(plot.yscale)
        p.grid()
        for args, kwargs in plot.curves: p.plot(*args, **kwargs)
        p.axes().set_aspect(plot.aspect)
        if plot.legend:
            p.rc('legend', fontsize='small')
            p.legend(shadow=0, loc='best')
        
        if not os.path.isdir(plot.dir): os.mkdir(plot.dir)
        if plot.pgf:
            p.savefig(plot.dir+plot.name+'.pgf')
            print(plot.name+'.pgf')
        if plot.pdf:
            p.savefig(plot.dir+plot.name+'.pdf', bbox_inches='tight')
            print(plot.name+'.pdf')
        
        p.close()
        
    def __reset(self):
        self.__plots = []
        self.__nr_show = 0
        self.__nr_id = 0