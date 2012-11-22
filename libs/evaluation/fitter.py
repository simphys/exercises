#!/usr/bin/python2
# -*- coding:utf-8 -*-

from __future__ import division
from scipy import optimize
import numpy as np

class Fitter(object):
    def __init__(self):
        pass
    def loadFunction(self, f, p0, maxfev=None):
        self.__calcedParams = False
        self.__calcedData = False
        self.__f = f
        self.__p0 = p0
        self.__maxfev = maxfev
    def loadData(self, x, y, sigma=None, scale = 'linear'):
        self.__calcedParams = False
        self.__calcedData = False
        self.__x = x
        self.__y = y
        self.__s = sigma
        if scale == 'log': self.__log = True
        else: self.__log = False
    def __call__(self, x):
        self.__calcParams()
        return self.__f(x, *self.__p)
    @property
    def params(self):
        self.__calcParams()
        return self.__p
    @property
    def std(self):
        self.__calcParams()
        return self.__std
    @property
    def correlation(self):
        self.__calcParams()
        return self.__correl
    @property
    def x(self):
        self.__calcData()
        return self.__xout
    @property
    def y(self):
        self.__calcData()
        return self.__yout
    def __calcParams(self):
        if not self.__calcedParams:
            if self.__maxfev==None: self.__maxfev = 100*(len(self.__x)+1)
            p, covar = optimize.curve_fit(self.__f, self.__x, self.__y, self.__p0, self.__s, maxfev = self.__maxfev)
            n = len(covar)
            std = np.zeros(n)
            correl = covar.copy()
            for i in range(n):
                std[i] = np.sqrt(covar[i][i])
            for i in range(n):
                correl[i,:] /= std[i]
                correl[:,i] /= std[i]
            self.__p = p
            self.__std = std
            self.__correl = correl
            self.__calcedParams = True
    def __calcData(self):
        self.__calcParams()
        if not self.__calcedData:
            if self.__log: self.__xout = np.logspace(np.log10(self.__x[0]), np.log10(self.__x[-1]), num=100)
            else: self.__xout = np.linspace(self.__x[0], self.__x[-1], num=500)
            self.__yout = self.__f(self.__xout, *self.__p)
            self.__calcedData = True