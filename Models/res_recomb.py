#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 18:06:28 2018

@author: usyd
"""
from pathlib import Path
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
#
## Generating and saving redshifts.
#zpicks = np.random.uniform(low=0.0001, high=1088, size=(1000,))
#zpicks[-1] = 1089
#zpicks = np.sort(zpicks, axis=None)
#pickle.dump(zpicks, open(f'zpicks_{zpicks[-1]}.p', 'wb'))

# Extracting pre-made redshifts z=0 to z=1089.
try:
    with open('data/zpicks_1000_1089.p','rb') as rfp: zpicks = pickle.load(rfp)
except:
    print("zpicks_1000_1089.p didnt't open")

data_dic = {'zpicks':zpicks}

## LCDM mag and da.
#names = ['Mcorr', 'matter']
#values = np.array([-19.3, 0.3])
#mmag, mda = datasim.magn(names, values, data_dic, 'LCDM', plot_key=False)

## Adding noise to LCDM mag and da.
#mag_mu, mag_sd = 0.0, 0.2
#nmag = datasim.gnoise(mmag, mag_mu, mag_sd)

timed = False

if timed:
    import cProfile, pstats, io
    pr = cProfile.Profile()
    pr.enable()

models = 'LCDM', 'stepfall', 'waterfall'
models = 'exotic', 'LCDM'
noise_options = 0.001, 0.07#, 0.14
npoints_options = 1048, 10480#, 104800
yaxis_tick = 1, 2, 3 # creating non real y values for visual separation of da's. da's plotted are all at z=1089.

#noise = 0.14
#npoints = 10480

for npoints in npoints_options:
    for noise in noise_options:
        da_list = []
        ml_da_list = []
        for test_key in models:
            if test_key:
                file_path = f'results_error_vs_data/{test_key}/sigma{noise}_npoints{npoints}.p'
                my_file = Path(file_path)
                if my_file.is_file():
                    # results of emcee runs
                    with open(file_path,'rb') as rfp: propert, sampler = pickle.load(rfp)
                else:
                    print(f"Couldn't open {file_path}")
                # Collecting every skip-th set of parameters and da they produce.
                skip = 100
                da_distrib = []
                for i in range(0, len(sampler.flatchain), skip):
                    values = sampler.flatchain[i, :] # no parsing real names as not needed
                    names = np.zeros(len(values)).tolist()
                    mag, da = datasim.magn(names, values, data_dic,
                                           test_key, plot_key=False)
                    da_distrib.append(da[-1])
                    ml_i = np.argmax(sampler.flatlnprobability) # max like index
                    ml_values = sampler.flatchain[ml_i,:]
                    ml_mag, ml_da = datasim.magn(names, ml_values, data_dic,
                                                test_key)
                da_list.append(da_distrib)
                ml_da_list.append(ml_da[-1])

            c = 'palevioletred', 'seagreen', 'slateblue', 'saddlebrown'
            ec = 'lightpink', 'lightgreen', 'lightblue', 'chocolate'

        plt.figure()
        plt.title(f'Angular diameter distances at z = {zpicks[-1]}')
        plt.ylabel('$(H_0 /c) * D_A$')
        plt.xlabel('z')
        plt.grid(True)
        plt.ylim(0.002,0.005)
        for i in range(len(models)):
            da_distrib = da_list[i]
            z_array = np.ones(len(da_distrib))*1089
            plt.scatter(z_array, da_distrib, s=40, facecolors='none', edgecolors=c[i], label=models[i])
        plt.legend()

        plt.figure()
        plt.title(f'Histogram of angular diameter distances at z = {zpicks[-1]}'+
        f'\n noise on data = {noise}, dataset size = {npoints}')
        plt.xlabel('$(H_0 /c) * D_A$')
        for i in range(len(models)):
            da_distrib = da_list[i]
            plt.hist(da_distrib, bins = 50, normed=True, color=c[i], histtype='step', stacked=True, fill=False, label=models[i])
        plt.legend()

        plt.figure()
        plt.title(f'$\mu$ and $\sigma$ of the angular diameter distance distribution'+
        f'\n noise on data = {noise}, dataset size = {npoints}')
        plt.xlabel('$(H_0 /c) * D_A$')
        plt.ylabel(r'$z = 1089$')
        for i in range(len(models)):
            da_distrib = np.asarray(da_list[i])
            da_mean = np.mean(da_distrib)
            da_sd = np.std(da_distrib)
            if da_sd < da_mean/100:
                print(f'da_sd = {da_sd}, model = {test_key}, mean/100 = {da_mean/100}')
            plt.errorbar(da_mean, yaxis_tick[i], xerr=da_sd, fmt='-o', color=c[i], ecolor=ec[i], elinewidth=3, capsize=0, label=models[i])
        plt.ylim(-4,8)
        plt.yticks([])
        plt.legend()

        plt.figure()
        plt.title(f'Max like and $\sigma$ of the angular diameter distance distribution'+
        f'\n noise on data = {noise}, dataset size = {npoints}')
        plt.xlabel('$(H_0 /c) * D_A$')
        plt.ylabel(r'$z = 1089$')
        for i in range(len(models)):
            da_distrib = np.asarray(da_list[i])
            da_sd = np.std(da_distrib)
            plt.errorbar(ml_da_list[i], yaxis_tick[i], xerr=da_sd, fmt='-o', color=c[i], ecolor=ec[i], elinewidth=3, capsize=0, label=models[i])
        plt.ylim(-4,8)
        plt.yticks([])
        plt.legend()
        plt.show()

if timed:
    pr.disable()
    s = io.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print (s.getvalue())

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