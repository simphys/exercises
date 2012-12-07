#!/usr/bin/python2
# -*- coding:utf-8 -*-

from __future__ import division
#import matplotlib
#matplotlib.use('GtkAgg')
import matplotlib as p
import numpy as np
import os

class Plot(object):
    def __init__(self,title,name,xlabel,ylabel,xscale,yscale,aspect,show,save,pgf,directory):
        self.curves = []
        self.title = title
        self.name = name
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.xscale = xscale
        self.yscale = yscale
        self.aspect = aspect
        self.show = show
        self.save = save
        self.pgf = pgf
        self.dir = directory
    def addCurve(self,xdata,ydata,fmt,label):
        self.curves.append([xdata,ydata,fmt,label])

class Plotter(object):
    def __init__(self, show = True, save = False, pgf = False, name='', directory=''):
        if not name: name = 'plot'
        self.__global_show = show
        self.__global_save = save
        self.__global_pgf = pgf
        self.__global_dir = directory
        self.__global_name = name
        self.__reset()
    def new(self, title='', name='', xlabel='', ylabel='', xscale='linear', yscale='linear', aspect='auto', show='nab', save='nab', pgf='nab', directory=''):
        if show == 'nab': show = self.__global_show
        if save == 'nab': save = self.__global_save
        if pgf == 'nab': pgf = self.__global_pgf
        if not title: title = name
        if not directory: directory = self.__global_dir
        self.__plots.append(Plot(title,name,xlabel,ylabel,xscale,yscale,aspect,show,save,pgf,directory))
        if show: self.__nr_show += 1
        if save: self.__nr_save += 1
    def plot(self, xdata, ydata=None, fmt='-',label='' ):
        self.__plots[-1].addCurve(xdata, ydata, fmt, label)
    def make(self, ncols = 2, swindow = (17,10), sfile = (8,3)):
        if self.__nr_show == 1: 
            f, axarr = p.subplots(1, 1, figsize=swindow)
            axarr = np.array([axarr])
        elif self.__nr_show > 1:
            nrows = int(np.ceil(1.*self.__nr_show/ncols))
            f, axarr = p.subplots(nrows, ncols, figsize=swindow)
        nplots = 0
        for n,plot in enumerate(self.__plots):
            if plot.show:
                self.__show(n,plot,axarr.flatten()[nplots])
                nplots += 1
            if plot.save: 
                self.__save(n,plot,sfile)
        if self.__nr_show > 0:
            f.tight_layout()
            p.show()
        self.__reset()
    def __show(self,n,plot, ax):
        ax.set_title("Plot %i) "%(n)+plot.title,fontsize='medium', fontweight='bold', x=.05, y =1., ha = 'left', va='bottom')
        ax.set_xlabel(plot.xlabel,fontsize='small')
        ax.set_ylabel(plot.ylabel,fontsize='small')
        ax.tick_params(labelsize='small')
        ax.set_xscale(plot.xscale)
        ax.set_yscale(plot.yscale)
        ax.grid()
        for curve in plot.curves:
            if curve[1] == None: ax.plot(curve[0], curve[2],label=curve[3])
            else: ax.plot(curve[0], curve[1], curve[2],label=curve[3])
        ax.legend(shadow=0, loc='best',fontsize='small')
        ax.set_aspect(plot.aspect)
    def __save(self,n,plot,sfile):
        p.figure(figsize=sfile)
        p.xlabel(plot.xlabel)
        p.ylabel(plot.ylabel)
        p.xscale(plot.xscale)
        p.yscale(plot.yscale)
        p.grid()
        for curve in plot.curves: 
            if curve[1] == None: p.plot(curve[0],curve[2], label=curve[3])
            else: p.plot(curve[0],  curve[1], curve[2], label=curve[3])
        p.rc('legend', fontsize='small')
        p.legend(shadow=0, loc='best')
        p.axes().set_aspect(plot.aspect)
        if not plot.dir: plot.dir = './plots/'
        if not plot.name: plot.name = self.__global_name+'_%0*i'%(2,n)
        if not os.path.isdir(plot.dir): os.mkdir(plot.dir)
        if plot.pgf: p.savefig(plot.dir+plot.name+'.pgf')
        else: p.savefig(plot.dir+plot.name+'.pdf', bbox_inches='tight')
        p.close()
    def __reset(self):
        self.__plots = []
        self.__nr_show = 0
        self.__nr_save = 0
