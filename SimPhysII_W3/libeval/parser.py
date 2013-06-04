#!/usr/bin/python2
# -*- coding:utf-8 -*-

from __future__ import division
import numpy as np
import os

class Parser(object):
    def __init__(self, path = "./", extension = ".txt", use=[0,1], sortbyname = True):
        self.__date = []
        self.__name = []
        self.__data = []
        for root, _dirs, filenames in os.walk(path):
            for name in filenames:
                if os.path.splitext(name)[-1] == extension:
                    filepath = os.path.join(root, name)
                    with open(filepath, 'r') as datei:
                        self.__date.append(os.path.getctime(os.path.join(root,name)))
                        self.__name.append(os.path.splitext(name)[0].replace('_','_')) # replace('_','\_') is important for LaTEX
                        self.__data.append(np.array([]))
                        for zeile in datei:
                            try:
                                liste = np.array(zeile.replace(',',' ').split()).astype(float)[use]
                            except:
                                pass
                            else:
                                self.__data[-1] = np.append(self.__data[-1],liste)
                        self.__data[-1] = self.__data[-1].reshape((-1,len(use))).T
                        sort = np.argsort(self.__data[-1][0])
                        self.__data[-1][:] = self.__data[-1][:,sort]                      
        if sortbyname: sort = np.argsort(self.__name)
        else: sort = np.argsort(self.__date)
        self.__date = np.array(self.__date)[sort]
        self.__name = np.char.lstrip(np.array(self.__name)[sort],'0')
        self.__data = np.array(self.__data)[sort]
    def __getitem__(self, index):
        return self.__name[index], self.__data[index]
    @property
    def len(self):
        return len(self.__name)
    @property
    def names(self):
        return self.__name
    @names.setter
    def names(self,names):
        self.__name = names
    @property
    def datas(self):
        return self.__data
    @datas.setter
    def datas(self,datas):
        self.__data = datas