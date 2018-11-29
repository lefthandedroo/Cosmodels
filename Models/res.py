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

# Pantheon data:
import pandas as pd
pantheon = pd.read_csv('./data/lcparam_full_long.txt', sep=" ")

# Reading each txt file column of interest as numpy.ndarray
mag = pantheon.mb.values
zpicks = pantheon.zhel.values

# Stacking arrays together and sorting by accending redshift.
data = np.stack((mag,zpicks), axis=0)
data.sort(axis=-1)

mag_p = data[0]
zpicks = data[-1]
data_dic = {'zpicks':zpicks}

## Pantheon data plot.
#plt.figure()
#plt.title('Pantheon')
#plt.ylabel('Mag')
#plt.xlabel('redshift')
#plt.scatter(zpicks, mag)

test_key = 'stepfall'

## Generating redshifts.
#zpicks = np.random.uniform(low=0.0, high=2.5, size=(1000,))
#zpicks = np.sort(zpicks, axis=None)
#data_dic = {'zpicks':zpicks}

# LCDM mag and da.
names = ['Mcorr', 'matter']
values = np.array([-19.3, 0.3])
mag0, da0 = datasim.magn(names, values, data_dic, 'LCDM', plot_key=False)

# Generating noisy LCDM mag and da.
mag_mu, mag_sd = 0.0, 0.07
mag = datasim.gnoise(mag0, mag_mu, mag_sd)
da_mu, da_sd = 0.0, 0.001
da = datasim.gnoise(da0, da_mu, da_sd)

## Does test_key reduce to LCDM?
#names = ['Mcorr','matter', 'radiation', 'a_ombar', 'b_ombar', 'c_ombar',
#         'v_in', 'w_in', 'x_in', 'y_in', 'z_in']
#values = np.array([-19.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
#mag1, da1 = datasim.magn(names, values, data_dic, test_key, plot_key=False)


# Model for test_key with fluids but no interaction.
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


plt.figure()
plt.title('Generated SN Ia magnitudes '+'\n Noise parameters: $\mu =$ %s, $\sigma =$ %s'%(mag_mu, mag_sd))
plt.xlabel('z')
plt.scatter(zpicks, mag_p, label='pantheon SN Ia data', marker=',', s=1)
plt.scatter(zpicks, mag, label='noisy LCDM data', marker=',', s=1)
plt.plot(zpicks, mag0, label='LCDM model', color='red')
#plt.plot(zpicks, mag1, label='test_key in LCDM mode')
plt.legend()

# Mag from parameters with highest likelihood.
bi = np.argmax(sampler.flatlnprobability)
values = sampler.flatchain[bi,:]
mag3, da3 = datasim.magn(names, values, data_dic, test_key, plot_key=False)

# Mag from mean params from emcee parameter distributions.
values = np.zeros((len(names)))
for i in range(len(names)):
    values[i] = np.mean(sampler.flatchain[:,i])
#mag_mean, da_mean = datasim.magn(names, values, data_dic, test_key, plot_key=False)


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
values = flat_sorted[2,:]
mag4, da4 = datasim.magn(names, values, data_dic, test_key, plot_key=False)

# Mag from parameters with lowest likelihood.
worst=[]
for i in range(len(values)):
    worst.append(flat_sorted[-1,i])
values = np.asarray(worst)
mag5, da5 = datasim.magn(names, values, data_dic, test_key, plot_key=False)

plt.figure()
plt.title('SN Ia magnitudes')
plt.ylabel('Mag')
plt.xlabel('z')
plt.scatter(zpicks, mag, label='noisy LCDM data', marker=',', s=1)
plt.scatter(zpicks, mag_p, label='pantheon SN Ia data', marker=',', s=1)
#plt.scatter(zpicks, mag_mean, label='mean emcee parameters', marker=',', s=1)
#plt.plot(zpicks, mag0, label='LCDM')
#plt.plot(zpicks, mag2, label= test_key+' no interaction')
plt.plot(zpicks, mag3, label= test_key+' max likelihood')
#plt.plot(zpicks, mag4, label= test_key+' 2nd highest likelihood')
#plt.plot(zpicks, mag5, label= test_key+' lowest likelihood')
plt.legend()

m_p_diff = mag_p - mag0
mdata_diff = mag - mag0
#mmean_diff = mag - mag_mean
mbest_diff = mag0 - mag3
msecond_best_diff = mag0 - mag4
mworst_diff = mag0 - mag5

plt.figure()
plt.title('SN Ia magnitude residuals')
plt.ylabel('Mag')
plt.xlabel('z')
plt.scatter(zpicks, m_p_diff, label='pantheon SN Ia data-LCDM', marker=',', s=1)
plt.scatter(zpicks, mdata_diff, label='fake LCDM data-LCDM', marker=',', s=1)
#plt.scatter(zpicks, mmean_diff, label='LCDM - mean emcee parameters', marker=',', s=1)
plt.plot(zpicks, mbest_diff, label='LCDM - max likelihood')
#plt.plot(zpicks, msecond_best_diff, label='LCDM - 2nd highest likelihood')
#plt.plot(zpicks, mworst_diff, label='LCDM - lowest likelihood')
plt.legend()


plt.figure()
plt.title('Generated angular diameter distances'+'\n Noise parameters: $\mu =$ %s, $\sigma =$ %s'%(da_mu, da_sd))
plt.ylabel('$(H_0 /c) * D_A$')
plt.xlabel('z')
plt.scatter(zpicks, da, label='noisy LCDM data', marker=',', s=1)
plt.plot(zpicks, da0, label='LCDM model', color='red')
#plt.scatter(zpicks, da_mean, label='mean emcee parameters', marker=',', s=1)
#plt.plot(zpicks, da1, label='test_key in LCDM mode')
plt.legend()

ddata_diff = da - da0
#dmean_diff = da - da_mean
dbest_diff = da0 - da3
dsecond_best_diff = da0 - da4
dworst_diff = da0 - da5

plt.figure()
plt.title('Angular diameter distances')
plt.ylabel('$(H_0 /c) * D_A$')
plt.xlabel('z')
plt.scatter(zpicks, da, label='noisy LCDM data', marker=',', s=1)
#plt.plot(zpicks, da0, label='LCDM')
#plt.plot(zpicks, da2, label= test_key+ ' no interaction')
plt.plot(zpicks, da3, label= test_key+ ' max likelihood')
#plt.plot(zpicks, da4, label= test_key+ ' 2nd highest likelihood')
#plt.plot(zpicks, da5, label= test_key+ ' lowest likelihood')
plt.legend()

plt.figure()
plt.title('Angular diameter distance residuals')
plt.xlabel('z')
plt.ylabel('$(H_0 /c) * D_A$')
plt.scatter(zpicks, ddata_diff, label='fake LCDM data-LCDM', marker=',', s=1)
#plt.scatter(zpicks, dmean_diff, label='LCDM - mean emcee parameters', marker=',', s=1)
plt.plot(zpicks, dbest_diff, label='LCDM - max likelihood')
#plt.plot(zpicks, dsecond_best_diff, label='LCDM - 2nd highest likelihood')
#plt.plot(zpicks, dworst_diff, label='LCDM - lowest likelihood')
plt.legend()

plt.show()