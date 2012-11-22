#!/usr/bin/python2
# -*- coding:utf-8 -*-

from __future__ import division
import matplotlib
matplotlib.use('GtkAgg')
import matplotlib.pyplot as p
import numpy as np
import os

class Plot(object):
    def __init__(self,title,name,xlabel,ylabel,xscale,yscale,show,save,pgf,directory):
        self.curves = []
        self.title = title
        self.name = name
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.xscale = xscale
        self.yscale = yscale
        self.show = show
        self.save = save
        self.pgf = pgf
        self.dir = directory
    def addCurve(self,xdata,ydata,fmt,label):
        self.curves.append([xdata,ydata,fmt,label])

class Plotter(object):
    def __init__(self, show = True, save = False, pgf = False, directory = ''):
        self.__global_show = show
        self.__global_save = save
        self.__global_pgf = pgf
        self.__global_dir = directory
        self.__reset()
    def new(self, title='', name='', xlabel='x', ylabel='y', xscale='linear', yscale='linear', show='nan', save='nan', pgf='nan', directory='nan'):
        if title == 'nan': title = name
        if show == 'nan': show = self.__global_show
        if save == 'nan': save = self.__global_save
        if pgf == 'nan': pgf = self.__global_pgf
        if directory == 'nan': directory = self.__global_dir
        self.__plots.append(Plot(title,name,xlabel,ylabel,xscale,yscale,show,save,pgf,directory))
        if show: self.__nr_show += 1
        if save: self.__nr_save += 1
    def plot(self, xdata, ydata, fmt='+-',label='' ):
        self.__plots[-1].addCurve(xdata, ydata, fmt, label)
    def make(self, ncols = 2):
        if self.__nr_show == 1: 
            f, axarr = p.subplots(1, 1, figsize=(15,10))
            axarr = np.array([axarr])
        elif self.__nr_show > 1:
            nrows = int(np.ceil(1.*self.__nr_show/ncols))
            f, axarr = p.subplots(nrows, ncols, figsize=(15,10))
        nplots = 0
        for n,plot in enumerate(self.__plots):
            if plot.show:
                self.__show(n,plot,axarr.flatten()[nplots])
                nplots += 1
            if plot.save: 
                self.__save(n,plot)
        if self.__nr_show > 0:
            f.tight_layout()
            p.show()
        self.__reset()
    def __show(self,n,plot, ax):
        ax.set_title("Plot %i - "%(n)+plot.title,fontsize='medium', fontweight='bold', x=.05, y =1., ha = 'left', va='bottom')
        ax.set_xlabel(plot.xlabel,fontsize='small')
        ax.set_ylabel(plot.ylabel,fontsize='small')
        ax.tick_params(labelsize='small')
        ax.set_xscale(plot.xscale)
        ax.set_yscale(plot.yscale)
        ax.grid()
        for curve in plot.curves: ax.plot(curve[0], curve[1], curve[2],label=curve[3])
        ax.legend(shadow=0, loc='best',fontsize='small')
    def __save(self,n,plot):
        p.figure(figsize=(8,3))
        p.xlabel(plot.xlabel)
        p.ylabel(plot.ylabel)
        p.xscale(plot.xscale)
        p.yscale(plot.yscale)
        p.grid()
        for curve in plot.curves: p.plot(curve[0], curve[1], curve[2], label=curve[3])
        p.rc('legend', fontsize='small')
        p.legend(shadow=0, loc='best')
        if not plot.dir: plot.dir = './plots/'
        if not plot.name: plot.name = os.path.splitext(os.path.basename(__file__))[0]+'_%0*i'%(2,n)
        if not os.path.isdir(plot.dir): os.mkdir(plot.dir)
        if plot.pgf: p.savefig(plot.dir+plot.name+'.pgf')
        else: p.savefig(plot.dir+plot.name+'.pdf', bbox_inches='tight')
        p.close()
    def __reset(self):
        self.__plots = []
        self.__nr_show = 0
        self.__nr_save = 0