#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:50:29 2017

@author: BallBlueMeercat
"""
import numpy as np

#test_params = {'m':0.3, 'gamma':0, 'M':-19}
#
#ndim = len(test_params)
#nwalkers = int(ndim * 2)
#
#poslist = list(test_params.values())
#print('poslist =',poslist)
#pos = []
#for i in poslist:
#    pos.append(i)
#startpos = np.array(pos)
#print('startpost =',startpos)
#pos = [startpos + 0.01*np.random.randn(ndim) for i in range(nwalkers)]
#print('pos =',pos)

test_params = {'m':0.3, 'gamma':0, 'M':-19}

ndim = len(test_params)
nwalkers = int(ndim * 2)

poslist = list(test_params.values())
print('poslist =',poslist)
startpos = np.array(poslist)
print('startpost =',startpos)
pos = [startpos + 0.01*np.random.randn(ndim) for i in range(nwalkers)]
print('pos =',pos)