#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 13:14:13 2018

@author: BallBlueMeercat
"""
import matplotlib.pyplot as plt
import numpy as np
import datasim

#data_name = 'synth'
data_name = 'pantheon'

if data_name == 'pantheon':
    print('-----Using pantheon')
    # Pantheon data:
    import pandas as pd
    pantheon = pd.read_csv('./data/lcparam_full_long.txt', sep=" ")
    # Stacking arrays together and sorting by accending redshift.
    data = np.stack((pantheon.mb.values,pantheon.zhel.values), axis=0)
    data.sort(axis=-1)
    mag = data[0]
    zpicks = data[-1]
    data_dic = {'mag':mag, 'zpicks':zpicks}
elif data_name == 'synth':
    print('-----Using artificial data')
    # Loading artificial LCDM SN Ia data:
    from pathlib import Path
    import pickle
    dataname = f'data/1048_3.0_sigma_0.07.p'
    my_file = Path(dataname)
    if my_file.is_file():
        with open(dataname,'rb') as rfp: zpicks, mag = pickle.load(rfp)
    data_dic = {'mag':mag, 'zpicks':zpicks}

scolnic_Mcorr_zpicks = [
0.30012
,0.14761
,0.22853
,0.1388
,0.14343
,0.09201
,0.1568
,0.23851
,0.14949
,0.50718
,0.25249
,0.28397
,0.2368
,0.34092
,0.26388]


scolnic_Mcorr = [
21.71
,19.89
,20.97
,19.82
,19.88
,18.95
,20.12
,21.20
,20.15
,22.99
,21.21
,21.32
,21.02
,21.89
,21.06]

scolnic_data_mag = np.zeros(len(scolnic_Mcorr))

i=0
for z in scolnic_Mcorr_zpicks:
    index = np.argmax(zpicks == z)
    if index == 0:
        print(f"Didn't find an exact z = {z} match")
        index = np.argmax(zpicks > z)
        print(f'but found closest z = {zpicks[index]}')
#    print(index)
    scolnic_data_mag[i] = mag[index]
    i+=1

plt.figure()
plt.title(f'SN Ia mag from Scolnic and mag at Scolnic paper redshifts')
plt.scatter(scolnic_Mcorr_zpicks, scolnic_Mcorr, label='scolnic paper mag')
plt.scatter(scolnic_Mcorr_zpicks, scolnic_data_mag, label=f'{data_name} data at scolnic z')
plt.legend()



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


test_key0 = 'rLCDM'
sample0 = np.loadtxt(f'./results_Bfactor/{data_name}/2_model_{test_key0}/sample.txt')
sample_info0 = np.loadtxt(f'./results_Bfactor/{data_name}/2_model_{test_key0}/sample_info.txt')
transposed_sample_info0 = sample_info0.transpose()
maxlike_index0 = np.argmax(transposed_sample_info0[1]) # The 1st (staring column 0th) in sample_info is the log-likelihood value.
names0 = ['Mcorr', 'm_ombar']
values0 = sample0[maxlike_index0,:]
#print('values0',values0)
mag0, da0 = datasim.magn(names0, values0, data_dic, test_key0, plot_key=False)

test_key1 = 'stepfall'
sample1 = np.loadtxt(f'./results_Bfactor/{data_name}/2_model_{test_key1}/sample.txt')
sample_info1 = np.loadtxt(f'./results_Bfactor/{data_name}/2_model_{test_key1}/sample_info.txt')
transposed_sample_info1 = sample_info1.transpose()
maxlike_index1 = np.argmax(transposed_sample_info1[1])
names1 = ['Mcorr', 'matter', 'radiation', 'a_ombar', 'v_in', 'w_in', 'x_in']
values1 = sample1[maxlike_index1,:]
#print('values1',values1)
mag1, da1 = datasim.magn(names1, values1, data_dic, test_key1, plot_key=False)

scolnic_omega_m = [
0.304
,0.317
,0.289
,0.281
,0.301
,0.307]

for omega_m in scolnic_omega_m:
    test_key00 = 'LCDM'
    names00 = ['Mcorr', 'm_ombar']
    values00 = [-19.3, omega_m]
    #print('values00',values00)
    mag00, da00 = datasim.magn(names00, values00, data_dic, test_key00, plot_key=False)

    ## SN Ia plots:
    #plt.figure()
    #plt.title(f'SN Ia magnitudes, data = {data_name} '+'\n Noise parameters: $\mu = 0.0$, $\sigma = 0.07$')
    #plt.xlabel('z')
    #plt.ylabel('Mag')
    #plt.scatter(zpicks, data_dic['mag'], label=data_name, marker=',', s=1)
    #plt.plot(zpicks, mag0, label=test_key0)
    #plt.plot(zpicks, mag1, label=test_key1)
    #plt.plot(zpicks, mag00, label=test_key00)
    #plt.legend()

    # Residuals:
    plt.figure()
    plt.title(f'SN Ia magnitude residuals, data = {data_name}')
    plt.ylabel('Mag')
    plt.xlabel('z')
    #plt.xlim(0.5,1)
    plt.scatter(zpicks, data_dic['mag'] - mag00, label=f'{data_name} SN Ia data - scolnic omega_m {test_key00}', marker=',', s=1)
    plt.scatter(zpicks, mag0 - mag1, label=f'DNest {test_key0} - DNest {test_key1}', marker=',', s=1)
    plt.scatter(zpicks, mag00 - mag1, label='LCDM w scolnic omega_m - DNest {test_key1}', marker=',', s=1)
    plt.grid(True)
    plt.legend()


plt.show()