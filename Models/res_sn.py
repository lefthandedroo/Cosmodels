#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 13:14:13 2018

@author: BallBlueMeercat
"""
import matplotlib.pyplot as plt
import numpy as np
import datasim

# Pantheon data:
import pandas as pd
pantheon = pd.read_csv('./data/lcparam_full_long.txt', sep=" ")
# Stacking arrays together and sorting by accending redshift.
data = np.stack((pantheon.mb.values,pantheon.zhel.values), axis=0)
data.sort(axis=-1)
pmag = data[0]
zpicks = data[-1]
data_dic = {'zpicks':zpicks}


# LCDM model magnitude and da
mmag, mda = datasim.magn(['Mcorr', 'matter'], np.array([-19.3, 0.3]), data_dic, 'LCDM')

## Does test_key reduce to LCDM?
#names = ['Mcorr', 'matter', 'radiation', 'y_in', 'z_in']
#values = np.array([-19.3, 0.3, 0.0, 0.0, 0.0])
#mag0, da0 = datasim.magn(names, values, data_dic, test_key, plot_key=False)
#
#plt.figure()
#plt.plot(zpicks, mmag, label='LCDM')
#plt.plot(zpicks, mag0, label=f'{test_key} in LCDM mode')
#plt.plot(zpicks, mmag-mag0, label='residual')
#plt.legend()
#plt.show()

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
#mag0, da0 = datasim.magn(names, values, data_dic, test_key, plot_key=False)



test_key = 'stepfall'
sample = np.loadtxt('./results_Bfactor/no05/2_model_stepfall/sample.txt')
sample_info = np.loadtxt('./results_Bfactor/no05/2_model_stepfall/sample_info.txt')
transposed_sample_info = sample_info.transpose()
maxlike_index = np.argmax(transposed_sample_info[1])

names = ['Mcorr', 'matter', 'radiation', 'a_ombar', 'v_in', 'w_in', 'x_in']
values = sample[maxlike_index,:]
mag1, da1 = datasim.magn(names, values, data_dic, test_key, plot_key=False)


test_key = 'LCDM'
sample = np.loadtxt('./results_Bfactor/no05/2_model_LCDM/sample.txt')
sample_info = np.loadtxt('./results_Bfactor/no05/2_model_LCDM/sample_info.txt')
transposed_sample_info = sample_info.transpose()
maxlike_index = np.argmax(transposed_sample_info[1])
names = ['Mcorr', 'm_ombar']
values = sample[maxlike_index,:]
mag0, da0 = datasim.magn(names, values, data_dic, test_key, plot_key=False)


# SN Ia plots:
plt.figure()
plt.title('SN Ia magnitudes '+'\n Noise parameters: $\mu = 0.0$, $\sigma = 0.07$')
plt.xlabel('z')
plt.ylabel('Mag')
plt.scatter(zpicks, pmag, label='pantheon', marker=',', s=1)
plt.plot(zpicks, mag0, label='LCDM')
plt.plot(zpicks, mag1, label='stepfall')
plt.legend()
plt.show()


m_p_diff = pmag - mmag                          # pantheon - LCDM
mbest_stepfall_diff = mmag - mag1               # LCDM - stepfall
mbest_LCDM_diff = mmag - mag0                   # LCDM - DNest LCDM params

# Residuals:
plt.figure()
plt.title('SN Ia magnitude residuals')
plt.ylabel('Mag')
plt.xlabel('z')
plt.scatter(zpicks, m_p_diff, label='pantheon SN Ia data - LCDM', marker=',', s=1)
plt.scatter(zpicks, mbest_stepfall_diff, label='LCDM - stepfall', marker=',', s=1)
plt.scatter(zpicks, mbest_LCDM_diff, label='LCDM - DNest LCDM', marker=',', s=1)
plt.grid(True)
plt.legend()





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