#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:50:29 2017

@author: BallBlueMeercat
"""
import numpy.random as rng
import pickle
import matplotlib.pyplot as plt
import numpy as np


#a = np.arange(5)
#hist, bin_edges = np.histogram(a, density=True)

#def wrap(x, a, b):
#    assert b > a
#    return (x - a)%(b - a) + a
#
#
#N = 3
#
#ombar_wm = 0.3    # starting from 0.3
#ombar_wr = 0.025    # starting from 0.025
#
#density_left = 1.0 - ombar_wm - ombar_wr       # e_total/e_crit = 1
#
#lower = 0.0
#upper = density_left/ N
#
#ombar_w01 = rng.rand()      # random number between 0 and 1
#ombar_w01 = wrap(ombar_w01, lower, upper)
#ombar_w05 = rng.rand()
#ombar_w05 = wrap(ombar_w05, lower, upper)
#ombar_w08 = rng.rand()
#ombar_w08 = wrap(ombar_w08, lower, upper)
#ombar_wde = density_left - ombar_w01 - ombar_w05 - ombar_w08
#
##print('ombar_m = ',ombar_wm)
##print('ombar_r = ',ombar_wr)
##print('ombar_w(-0.1) = ',ombar_w01)
##print('ombar_w(-0.5) = ',ombar_w05)
##print('ombar_w(-0.8) = ',ombar_w08)
##print('ombar_de = ',ombar_wde)
#
#total = ombar_wm +ombar_wr +ombar_w01 +ombar_w05 +ombar_w08 +ombar_wde
#print('total =',total)
#
#import datasim
#import numpy as np
#try:
#    with open('zpicks_1089.p','rb') as rfp: zpicks = pickle.load(rfp)
#except:
#    print("zpicks_1089.p didnt't open")
#print(len(zpicks))
#data_dic = {'zpicks':zpicks}
#mag, da = datasim.magn(['Mcorr', 'matter'], np.array([-19.3, 0.3]), data_dic, 'LCDM', plot_key=False)
#plt.figure()
#plt.title(f'Artificial data, zpicks_1089')
#plt.scatter(zpicks, mag)
#


try:
    with open('data/zpicks_1048_3.p','rb') as rfp: all_zpicks = pickle.load(rfp)
except:
    print("zpicks_1048_3.p didnt't open")
#
#data_dic = {'zpicks':zpicks}
#mag, da = datasim.magn(['Mcorr', 'matter'], np.array([-19.3, 0.3]), data_dic, 'LCDM', plot_key=False)
#plt.figure()
#plt.title(f'Artificial data, zpicks_10890')
#plt.scatter(zpicks, mag)
#
#plt.show()
#
## Generating and saving redshifts.
#zpicks = np.random.uniform(low=0.0001, high=3, size=(1048,))
#zpicks = np.sort(zpicks, axis=None)
#zpicks[-1] = 3
#pickle.dump(zpicks, open(f'zpicks_{len(zpicks)}_{zpicks[-1]}.p', 'wb'))