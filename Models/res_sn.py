#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 13:14:13 2018

@author: BallBlueMeercat
"""
import matplotlib.pyplot as plt
import numpy as np
import datasim

data_name = 'synth'
data_name = 'pantheon'

if data_name == 'pantheon':
    # Pantheon data:
    import pandas as pd
    pantheon = pd.read_csv('./data/lcparam_full_long.txt', sep=" ")
    # Stacking arrays together and sorting by accending redshift.
    data = np.stack((pantheon.mb.values,pantheon.zhel.values), axis=0)
    data.sort(axis=-1)
    pmag = data[0]
    zpicks = data[-1]
    data_dic = {'mag':pmag, 'zpicks':zpicks}
elif data_name == 'synth':
    # Loading artificial LCDM SN Ia data:
    from pathlib import Path
    import pickle
    dataname = f'data/1048_3.0_sigma_0.07.p'
    my_file = Path(dataname)
    if my_file.is_file():
        with open(dataname,'rb') as rfp: zpicks, mag = pickle.load(rfp)
    data_dic = {'mag':mag, 'zpicks':zpicks}


# LCDM model magnitude and da
#mmag, mda = datasim.magn(['Mcorr', 'matter'], np.array([-19.3, 0.3]), data_dic, 'LCDM')

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
#mag, da = datasim.magn(names, values, data_dic, test_key, plot_key=False)


test_key1 = 'stepfall'
sample = np.loadtxt(f'./results_Bfactor/{data_name}/2_model_{test_key1}/sample.txt')
sample_info = np.loadtxt(f'./results_Bfactor/{data_name}/2_model_{test_key1}/sample_info.txt')
transposed_sample_info = sample_info.transpose()
maxlike_index = np.argmax(transposed_sample_info[1])

names = ['Mcorr', 'matter', 'radiation', 'a_ombar', 'v_in', 'w_in', 'x_in']
values = sample[maxlike_index,:]
mag1, da1 = datasim.magn(names, values, data_dic, test_key1, plot_key=False)


test_key0 = 'LCDM'
test_key0 = 'rLCDM'
sample = np.loadtxt(f'./results_Bfactor/{data_name}/2_model_{test_key0}/sample.txt')
sample_info = np.loadtxt(f'./results_Bfactor/{data_name}/2_model_{test_key0}/sample_info.txt')
transposed_sample_info = sample_info.transpose()
maxlike_index = np.argmax(transposed_sample_info[1])

names = ['Mcorr', 'm_ombar']
values = sample[maxlike_index,:]
mag0, da0 = datasim.magn(names, values, data_dic, test_key0, plot_key=False)

# SN Ia plots:
plt.figure()
plt.title(f'SN Ia magnitudes, data = {data_name} '+'\n Noise parameters: $\mu = 0.0$, $\sigma = 0.07$')
plt.xlabel('z')
plt.ylabel('Mag')
plt.scatter(zpicks, data_dic['mag'], label=data_name, marker=',', s=1)
plt.plot(zpicks, mag0, label=test_key0)
plt.plot(zpicks, mag1, label=test_key1)
plt.legend()

# Residuals:
plt.figure()
plt.title(f'SN Ia magnitude residuals, data = {data_name}')
plt.ylabel('Mag')
plt.xlabel('z')
#plt.xlim(0.5,1)
plt.scatter(zpicks, data_dic['mag'] - mag0, label=f'{data_name} SN Ia data - DNest {test_key0}', marker=',', s=1)
plt.scatter(zpicks, mag0 - mag1, label=f'DNest {test_key0} - DNest {test_key1}', marker=',', s=1)
#plt.scatter(zpicks, mbest_LCDM_diff, label='LCDM - DNest {test_key0}', marker=',', s=1)
plt.grid(True)
plt.legend()


plt.show()