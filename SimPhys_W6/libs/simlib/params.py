#!/usr/bin/python2
# -*- coding:utf-8 -*-
# Author: Sebastian Weber

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import md5

class Params(object):
    def __setattr__(self, name, value):
        object.__setattr__(self, name, value)
        if name != '_params__list': print("set %s\t= %s" % (name, repr(value)))
    def __getattr__(self, name):
        raise TypeError("Parameter '%s' is not set!" % (name))
    def id(self):
        return md5.new(repr(self.__dict__)).hexdigest()   
    def args(self):
        args = ''
        for k, v in self.__dict__.items(): args += ' --'+k+' %s'%(v)
        return args