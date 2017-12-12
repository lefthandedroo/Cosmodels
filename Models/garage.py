#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:50:29 2017

@author: BallBlueMeercat
"""

#import numpy as np
#import random
#from math import log10


start = 0
stop = 2
step = 0.05
i = 0
zlist = [start]


i = 0
#while True:
#    print('zlist is: ',zlist)
#    if zlist[i] > stop:
#        break
#    else:
#        nextvalue = zlist[i] + step
#        zlist.append(nextvalue)
#    i += 1


while zlist[i] < stop:
    nextvalue = zlist[i] + step
    zlist.append(nextvalue)
    i += 1
    continue
    print('i is ', i)
    print('zlist[i] is ',zlist[i])
print('zlist from zlist is: ',zlist)

#M= -19
#n = 10
#
#dlmpc = [0.86446572931294097,
# 0.64606397582287411,
# 0.0048361872398661361,
# 0.16798531977597295,
# 0.24523827944403684,
# 0.49173616689193961,
# 0.29433452566665375,
# 0.029537178192839499,
# 1.0414629022019111,
# 0.37053687119342366,
# 0.45]
#
#z = [5370995727.2544966,
# 3758068614.1052346,
# 20227789.392998893,
# 785777353.18180418,
# 1199628706.0813937,
# 2701836427.4300051,
# 1477914264.6854413,
# 125869569.09496653,
# 6760294905.2611141,
# 1931475800.7488716,
# 3]
#
##dlmpc = np.asarray(dlmpc)
##z = np.asarray(z)
##z.flatten()
#
##def msim(dlmpc, z, n):
##    """
##    Takes in:
##        dlmpc = luminosity distances in Mpc;
##        z = redshifts, dimensionless;
##        n = dimensionless number of data points to be generated.
##    Returns n apparent magnitudes m with redshit z < 2.
##    """
#zmax = 2        # Largest meausured z for a supernovae is 2.
#maxindex_i = np.where(z > zmax) 
#maxindex_i = np.asarray(maxindex_i)   # Converting to np array.
#   
#if maxindex_i.any():              # Check if instances of z > zmax exist.   
#    max_index = maxindex_i[0,0]
#else:
#    print('msim found no z above z = %s'%(zmax))
#
#dlmpc = dlmpc[:max_index]
#z = z[:max_index]
#  
#index_opts = range(len(z))
#data_ind = random.sample(index_opts, n)
#
#listdlmpc = []
#listz = []
#
#for i in data_ind:
#    listdlmpc.append(dlmpc[i])
#    listz.append(z[i])
#
#dlmpc = np.asarray(listdlmpc)
#z = np.asarray(listz)
#
## Calculatig apparent magnitudes of supernovae at the simulated
## luminosity distances using the distance modulus formula.
#m = []
#for i in range(len(dlmpc)):
#    m.append(5 * log10(dlmpc[i]/10) + M) 
#
##    return z, m
##
##z, m = msim(dlmpc, z, n)