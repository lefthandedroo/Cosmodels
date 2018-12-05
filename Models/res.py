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

## Pantheon data:
#import pandas as pd
#pantheon = pd.read_csv('./data/lcparam_full_long.txt', sep=" ")
## Stacking arrays together and sorting by accending redshift.
#data = np.stack((pantheon.mb.values,pantheon.zhel.values), axis=0)
#data.sort(axis=-1)
#pmag = data[0]
#zpicks = data[-1]
#data_dic = {'zpicks':zpicks}
## Pantheon data plot.
#plt.figure()
#plt.title('Pantheon')
#plt.ylabel('Mag')
#plt.xlabel('redshift')
#plt.scatter(zpicks, pmag, marker=',', s=1)

# Generating redshifts.
zpicks = np.random.uniform(low=0.0, high=1088, size=(1000,))
zpicks[-1] = 1089
zpicks = np.sort(zpicks, axis=None)
data_dic = {'zpicks':zpicks}

# LCDM mag and da.
names = ['Mcorr', 'matter']
values = np.array([-19.3, 0.3])
mmag, mda = datasim.magn(names, values, data_dic, 'LCDM', plot_key=False)

# Adding noise to LCDM mag and da.
mag_mu, mag_sd = 0.0, 0.2
nmag = datasim.gnoise(mmag, mag_mu, mag_sd)

test_key = 'exotic'

## Does test_key reduce to LCDM?
#names = ['Mcorr','matter', 'radiation', 'a_ombar', 'b_ombar', 'c_ombar',
#         'v_in', 'w_in', 'x_in', 'y_in', 'z_in']
#values = np.array([-19.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
#mag0, da0 = datasim.magn(names, values, data_dic, test_key, plot_key=False)

## Model for test_key with fluids but no interaction.
#if test_key == 'waterfall':
#    names = ['Mcorr',
#             'matter', 'radiation', 'a_ombar', 'b_ombar', 'c_ombar',
#             'v_in', 'w_in', 'x_in', 'y_in', 'z_in']
#    values = np.array([-19.3,
#                       0.3, 0.025, 0.1, 0.1, 0.1,
#                       0.0, 0.0, 0.0, 0.0, 0.0])
#elif test_key == 'stepfall':
#    names = ['Mcorr', 'matter', 'radiation', 'a_ombar',
#             'v_in', 'w_in', 'x_in']
#    values = np.array([-19.3, 0.3, 0.025, 0.1, 0.0, 0.0, 0.0])
#elif test_key == 'exotic':
#    names = ['Mcorr', 'matter', 'radiation', 'gamma', 'zeta']
#    values = np.array([-19.3, 0.3, 0.025, 0.0, 0.0])
#mag1, da1 = datasim.magn(names, values, data_dic, test_key, plot_key=False)

# Sampler from chosen test_key run.
# sd 0.07

from pathlib import Path

file_path = 'results_emcee/long/0.07_'+test_key+'/sampler.p'
my_file = Path(file_path)
if my_file.is_file():
    with open(file_path,'rb') as rfp:
        sampler007 = pickle.load(rfp)
else:
    file_path = 'results_emcee/long/0.07_'+test_key+'/sampler.p.zip'
    with open(file_path,'rb') as rfp:
        sampler007 = pickle.load(rfp)

flatlnprobability007 = sampler007.flatlnprobability
transposed_flatchain = sampler007.flatchain.transpose()
flat_sorted007 = np.vstack([transposed_flatchain, flatlnprobability007])
flat_sorted007.sort(axis=-1)

da_list = []
for i in range(len(sampler007.flatchain)):
    values = sampler007.flatchain[i, :]
    mag, da = datasim.magn(names, values, data_dic, test_key, plot_key=False)
    da_list.append(da)
    i += 10

plt.figure()
plt.title('Angular diameter distances')
plt.ylabel('$(H_0 /c) * D_A$')
plt.xlabel('z')
for da in da_list:
    plt.plot(zpicks, da)
plt.grid(True)
plt.legend()

## Mag from parameters with max likelihood.
#bi_007 = np.argmax(sampler007.flatlnprobability)
#values = sampler007.flatchain[bi_007,:]
#mag_max_007, da_max_007 = datasim.magn(names, values, data_dic, test_key, plot_key=False)
#
## Mag from parameters with 2nd highest likelihood.
#values = flat_sorted007[2,:]
#mag_2max_007, da_2max_007 = datasim.magn(names, values, data_dic, test_key, plot_key=False)
#
### Mag from parameters with lowest likelihood.
##values = flat_sorted007[-1,:]
##mag_min_007, da_min_007 = datasim.magn(names, values, data_dic, test_key, plot_key=False)

## sd 0.2
#with open('results_emcee/long/0.2_'+test_key+'/sampler.p','rb') as rfp:
#    sampler02 = pickle.load(rfp)
#
#flatlnprobability02 = sampler02.flatlnprobability
#transposed_flatchain = sampler02.flatchain.transpose()
#flat_sorted02 = np.vstack([transposed_flatchain, flatlnprobability02])
#flat_sorted02.sort(axis=-1)
#
#da_list = []
#for i in range(len(sampler02.flatchain)):
#    values = sampler02.flatchain[i, :]
#    mag, da = datasim.magn(names, values, data_dic, test_key, plot_key=False)
#    da_list.append(da)
#
#plt.figure()
#plt.title('Angular diameter distances')
#plt.ylabel('$(H_0 /c) * D_A$')
#plt.xlabel('z')
#for da in da_list:
#    plt.plot(zpicks, da, color='red')
#plt.grid(True)
#plt.legend()

## Mag from parameters with max likelihood.
#bi_02 = np.argmax(sampler02.flatlnprobability)
#values = sampler02.flatchain[bi_02,:]
#mag_max_02, da_max_02 = datasim.magn(names, values, data_dic, test_key, plot_key=False)
#
## Mag from parameters with 2nd highest likelihood.
#values = flat_sorted02[2,:]
#mag_2max_02, da_2max_02 = datasim.magn(names, values, data_dic, test_key, plot_key=False)
#
### Mag from parameters with lowest likelihood.
##values = flat_sorted02[-1,:]
##mag_min_02, da_min_02 = datasim.magn(names, values, data_dic, test_key, plot_key=False)
#
#
## SN Ia plots:
#plt.figure()
#plt.title('SN Ia magnitudes '+'\n Noise parameters: $\mu =$ %s, $\sigma =$ %s'%(mag_mu, mag_sd))
#plt.xlabel('z')
##plt.scatter(zpicks, pmag, label='pantheon', marker=',', s=1)
#plt.scatter(zpicks, nmag, label='noisy LCDM', marker=',', s=1)
#plt.plot(zpicks, mmag, label='LCDM', color='red')
##plt.plot(zpicks, mag0, label=test_key+' in LCDM mode')
##plt.plot(zpicks, mag1, label=test_key+' w no interaction')
#plt.plot(zpicks, mag_max_007, label= test_key+' max likelihood, $\sigma = 0.07$')
#plt.plot(zpicks, mag_max_02, label= test_key+' max likelihood, $\sigma = 0.2$')
##plt.plot(zpicks, mag_2max_007, label= test_key+' 2nd highest likelihood, $\sigma = 0.07$')
##plt.plot(zpicks, mag_2max_02, label= test_key+' 2nd highest likelihood, $\sigma = 0.2$')
##plt.plot(zpicks, mag_min_007, label= test_key+' lowest likelihood, $\sigma = 0.07$')
##plt.plot(zpicks, mag_min_02, label= test_key+' lowest likelihood, $\sigma = 0.2$')
#plt.grid(True)
#plt.legend()
#
##m_p_diff = pmag - mmag                          # pantheon - LCDM
#nmag_diff = nmag - mmag                         # noisy LCDM - LCDM
#mbest_diff_007 = mmag - mag_max_007             # LCDM - max like sd=0.07
#msecond_best_diff_007 = mmag - mag_2max_007     # LCDM - 2 highest like sd=0.07
##mworst_diff_007 = mmag - mag_min_007            # LCDM - lowest like sd=0.07
#mbest_diff_02 = mmag - mag_max_02               # LCDM - max like sd=0.2
#msecond_best_diff_02 = mmag - mag_2max_02       # LCDM - 2 highest like sd=0.2
##mworst_diff_02 = mmag - mag_min_02              # LCDM - lowest like sd=0.2
#
## Residuals:
#plt.figure()
#plt.title('SN Ia magnitude residuals')
#plt.ylabel('Mag')
#plt.xlabel('z')
##plt.scatter(zpicks, m_p_diff, label='pantheon SN Ia data-LCDM', marker=',', s=1)
#plt.scatter(zpicks, nmag_diff, label='noisy LCDM data-LCDM', marker=',', s=1)
#plt.plot(zpicks, mbest_diff_007, label='LCDM - max likelihood, $\sigma = 0.07$')
##plt.plot(zpicks, msecond_best_diff_007, label='LCDM - 2nd highest likelihood, $\sigma = 0.07$')
##plt.plot(zpicks, mworst_diff_007, label='LCDM - lowest likelihood, $\sigma = 0.07$')
#plt.plot(zpicks, mbest_diff_02, label='LCDM - max likelihood, $\sigma = 0.2$')
##plt.plot(zpicks, msecond_best_diff_02, label='LCDM - 2nd highest likelihood, $\sigma = 0.2$')
##plt.plot(zpicks, mworst_diff_02, label='LCDM - lowest likelihood, $\sigma = 0.2$')
#plt.grid(True)
#plt.legend()
#
#
#
#
#
## Angular diameter distance plots:
#plt.figure()
#plt.title('Angular diameter distances')
#plt.ylabel('$(H_0 /c) * D_A$')
#plt.xlabel('z')
#plt.plot(zpicks, mda, label='LCDM', color='red')
##plt.plot(zpicks, da0, label=test_key+' in LCDM mode')
##plt.plot(zpicks, da1, label=test_key+' w no interaction')
#plt.plot(zpicks, da_max_007, label= test_key+' max likelihood, $\sigma = 0.07$')
#plt.plot(zpicks, da_max_02, label= test_key+' max likelihood, $\sigma = 0.2$')
##plt.plot(zpicks, da_2max_007, label= test_key+' 2nd highest likelihood, $\sigma = 0.07$')
##plt.plot(zpicks, da_2max_02, label= test_key+' 2nd highest likelihood, $\sigma = 0.2$')
##plt.plot(zpicks, da_min_007, label= test_key+' lowest likelihood, $\sigma = 0.07$')
##plt.plot(zpicks, da_min_02, label= test_key+' lowest likelihood, $\sigma = 0.2$')
#plt.grid(True)
#plt.legend()
#
#dabest_diff_007 = mda - da_max_007             # LCDM - max like sd=0.07
#dasecond_best_diff_007 = mda - da_2max_007     # LCDM - 2 highest like sd=0.07
##daworst_diff_007 = mda - da_min_007            # LCDM - lowest like sd=0.07
#dabest_diff_02 = mda - da_max_02               # LCDM - max like sd=0.2
#dasecond_best_diff_02 = mda - da_2max_02       # LCDM - 2 highest like sd=0.2
##daworst_diff_02 = mda - da_min_02              # LCDM - lowest like sd=0.2
#
## Residuals:
#plt.figure()
#plt.title('Angular diameter distance residuals')
#plt.ylabel('$(H_0 /c) * D_A$')
#plt.xlabel('z')
#plt.plot(zpicks, dabest_diff_007, label='LCDM - max likelihood, $\sigma = 0.07$')
##plt.plot(zpicks, dasecond_best_diff_007, label='LCDM - 2nd highest likelihood, $\sigma = 0.07$')
##plt.plot(zpicks, daworst_diff_007, label='LCDM - lowest likelihood, $\sigma = 0.07$')
#plt.plot(zpicks, dabest_diff_02, label='LCDM - max likelihood, $\sigma = 0.2$')
##plt.plot(zpicks, dasecond_best_diff_02, label='LCDM - 2nd highest likelihood, $\sigma = 0.2$')
##plt.plot(zpicks, daworst_diff_02, label='LCDM - lowest likelihood, $\sigma = 0.2$')
#plt.grid(True)
#plt.legend()


plt.show()