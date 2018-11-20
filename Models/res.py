#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 18:06:28 2018

@author: usyd
"""
import pickle
import numpy as np
import matplotlib.pyplot as plt
import datasim


# Statistical parameteres of noise: mean, standard deviation.
mu, sigma = 0.0, 0.07 # sigma != 0

## Pantheon data:
#import pandas as pd
#pantheon = pd.read_csv('./data/lcparam_full_long.txt', sep=" ")
#
## Reading each txt file column of interest as numpy.ndarray
#mag = pantheon.mb.values
#zpicks = pantheon.zhel.values
#
## Stacking arrays together and sorting by accending redshift.
#data = np.stack((mag,zpicks), axis=0)
#data.sort(axis=-1)
#
#mag = data[0]
#zpicks = data[-1]
#data_dic = {'mag':mag, 'zpicks':zpicks}
#
## Pantheon data plot.
#plt.figure()
#plt.title('Pantheon')
#plt.ylabel('Mag')
#plt.xlabel('redshift')
#plt.scatter(zpicks, mag)

test_key = 'stepfall'

# Generating redshifts.
zpicks = np.random.uniform(low=0.0, high=100, size=(1000,))
zpicks = np.sort(zpicks, axis=None)
data_dic = {'mag':None, 'zpicks':zpicks}

# Generating noisy LCDM mag and da.
names = ['Mcorr', 'matter']
values = np.array([-19.3, 0.3])
mag, da = datasim.noisy_mag(mu, sigma, names, values, data_dic, 'LCDM')
data_dic = {'mag':mag, 'zpicks':zpicks}

# LCDM mag and da.
names = ['Mcorr', 'matter']
values = np.array([-19.3, 0.3])
mag0, da0 = datasim.magn(names, values, data_dic, 'LCDM', plot_key=False)

## Mag and da for test_key in LCDM mode.
#names = ['Mcorr','matter', 'radiation', 'a_ombar', 'b_ombar', 'c_ombar',
#         'v_in', 'w_in', 'x_in', 'y_in', 'z_in']
#values = np.array([-19.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
#mag1, da1 = datasim.magn(names, values, data_dic, test_key, plot_key=False)

plt.figure()
plt.title('SN Ia magnitudes')
plt.xlabel('z')
plt.scatter(zpicks, mag, label='noisy LCDM data', marker=',', s=1)
plt.plot(zpicks, mag0, label='LCDM model', color='orange')
#plt.plot(zpicks, mag1, label='test_key in LCDM mode')
plt.legend()

plt.figure()
plt.title('Angular diameter distances')
plt.ylabel('(H_0/c) * D_A')
plt.xlabel('z')
plt.scatter(zpicks, da, label='noisy LCDM data', marker=',', s=1)
plt.plot(zpicks, da0, label='LCDM model', color='orange')
#plt.plot(zpicks, da1, label='test_key in LCDM mode')
plt.legend()

# Mag for test_key with fluids but no interaction.
if test_key == 'waterfall':
    names = ['Mcorr',
             'matter', 'radiation', 'a_ombar', 'b_ombar', 'c_ombar',
             'v_in', 'w_in', 'x_in', 'y_in', 'z_in']
    values = np.array([-19.3,
                       0.3, 0.025, 0.1, 0.1, 0.1,
                       0.0, 0.0, 0.0, 0.0, 0.0])
elif test_key == 'stepfall':
    names = ['Mcorr', 'matter', 'radiation', 'a_ombar',
             'v_in', 'w_in', 'x_in']
    values = np.array([-19.3, 0.3, 0.025, 0.1, 0.0, 0.0, 0.0])
elif test_key == 'exotic':
    names = ['Mcorr', 'matter', 'radiation', 'gamma', 'zeta']
    values = np.array([-19.3, 0.3, 0.025, 0.0, 0.0])
mag2, da2 = datasim.magn(names, values, data_dic, test_key, plot_key=False)

# Sampler from chosen test_key run.
with open('results_emcee/long/'+test_key+'/sampler.p','rb') as rfp:
    sampler = pickle.load(rfp)

# Mag from parameters with highest likelihood.
bi = np.argmax(sampler.flatlnprobability)
values = sampler.flatchain[bi,:]
mag3, da3 = datasim.magn(names, values, data_dic, test_key, plot_key=False)

# Parameters with X best likelihood.
flatlnprobability = sampler.flatlnprobability
flatchain_M = sampler.flatchain[:,0]
flatchain_m = sampler.flatchain[:,1]
flatchain_r = sampler.flatchain[:,2]
flatchain_a = sampler.flatchain[:,3]
#flatchain_b = sampler.flatchain[:,4]
#flatchain_c = sampler.flatchain[:,5]
flatchain_v = sampler.flatchain[:,4]
flatchain_w = sampler.flatchain[:,5]
flatchain_x = sampler.flatchain[:,6]
#flatchain_y = sampler.flatchain[:,9]
#flatchain_z = sampler.flatchain[:,10]
# Stacking them together and sorting by accending redshift.
flat_sorted = np.stack((flatchain_M, flatchain_m, flatchain_r,
                        flatchain_a, flatchain_v, flatchain_w,
                        flatchain_x,flatlnprobability), axis=0)
flat_sorted.sort(axis=-1)

# Mag from parameters with second highest likelihood.
second_best=[]
for i in range(len(values)):
    second_best.append(sampler.flatchain[2,i])
values = np.asarray(second_best)
mag4, da4 = datasim.magn(names, values, data_dic, test_key, plot_key=False)

# Mag from parameters with lowest likelihood.
worst=[]
for i in range(len(values)):
    worst.append(sampler.flatchain[-2,i])
values = np.asarray(worst)
mag5, da5 = datasim.magn(names, values, data_dic, test_key, plot_key=False)

plt.figure()
plt.title('SN Ia magnitudes')
plt.ylabel('Mag')
plt.xlabel('z')
plt.scatter(zpicks, mag, label='noisy LCDM data', marker=',', s=1)
plt.plot(zpicks, mag0, label='LCDM')
plt.plot(zpicks, mag2, label= test_key+ ' no interaction')
plt.plot(zpicks, mag3, label= test_key+ ' max likelihood')
plt.plot(zpicks, mag4, label= test_key+ ' second highest likelihood')
plt.plot(zpicks, mag5, label= test_key+ ' lowest likelihood')
plt.legend()


mdata_diff = mag - mag0
mbest_diff = mag0 - mag3
msecond_best_diff = mag0 - mag4
mworst_diff = mag0 - mag5

plt.figure()
plt.title('Mag residuals')
plt.ylabel('Mag')
plt.xlabel('z')
plt.scatter(zpicks, mdata_diff, label='noisy LCDM data-LCDM', marker=',', s=1)
plt.plot(zpicks, mbest_diff, label='LCDM - best emcee fit')
plt.plot(zpicks, msecond_best_diff, label='LCDM - 2nd best emcee fit')
plt.plot(zpicks, mworst_diff, label='LCDM - worst emcee fit')
plt.legend()


ddata_diff = da - da0
dbest_diff = da0 - da3
dsecond_best_diff = da0 - da4
dworst_diff = da0 - da5

plt.figure()
plt.title('Angular diameter distances')
plt.ylabel('(H_0/c) * D_A')
plt.xlabel('z')
plt.scatter(zpicks, da, label='noisy LCDM data', marker=',', s=1)
plt.plot(zpicks, da0, label='LCDM')
plt.plot(zpicks, da2, label= test_key+ ' no interaction')
plt.plot(zpicks, da3, label= test_key+ ' max likelihood')
plt.plot(zpicks, da4, label= test_key+ ' second highest likelihood')
plt.plot(zpicks, da5, label= test_key+ ' lowest likelihood')
plt.legend()

plt.figure()
plt.title('DA residuals')
plt.xlabel('z')
plt.ylabel('(H_0/c) * D_A')
plt.scatter(zpicks, ddata_diff, label='noisy LCDM data-LCDM', marker=',', s=1)
plt.plot(zpicks, dbest_diff, label='LCDM - best emcee fit')
plt.plot(zpicks, dsecond_best_diff, label='LCDM - 2nd best emcee fit')
plt.plot(zpicks, dworst_diff, label='LCDM - worst emcee fit')
plt.legend()

plt.show()