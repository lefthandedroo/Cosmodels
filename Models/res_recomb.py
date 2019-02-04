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

zpicks = np.array([1089])
data_dic = {'zpicks':zpicks}

timed = False

if timed:
    import cProfile, pstats, io
    pr = cProfile.Profile()
    pr.enable()

models = 'exotic', 'LCDM'
noise_options = 0.01, 0.07, 0.14
npoints_options = 1048, 10480, 104800
# creating non-real y values for visual separation of da's. da's plotted are all at z=1089.
yaxis_tick = 1, 2, 3
#noise = 0.14
#npoints = 10480

for npoints in npoints_options:
    for noise in noise_options:
        da_list = []
        ml_da_list = [] # max likelihood da's
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
                    values = sampler.flatchain[i, :]
                    # not parsing real names as not needed
                    names = np.zeros(len(values)).tolist()
                    mag, da = datasim.magn(names, values, data_dic,
                                           test_key, plot_key=False)
                    da_distrib.append(da[-1]) # adding da at z[-1](=1089)
                    ml_i = np.argmax(sampler.flatlnprobability) # max like index
                    ml_values = sampler.flatchain[ml_i,:]
                    ml_mag, ml_da = datasim.magn(names, ml_values, data_dic,
                                                test_key)
                da_list.append(da_distrib)
                ml_da_list.append(ml_da[-1])

            c = 'palevioletred', 'seagreen', 'slateblue', 'saddlebrown'
            ec = 'lightpink', 'lightgreen', 'lightblue', 'chocolate'

        plt.figure()
        plt.title(f'Angular diameter distances at z = {zpicks[-1]},'
                +f'\n $\sigma$ on data={noise}, number of SN Ia used = {npoints}')
        plt.ylabel('$(H_0 /c) * D_A$')
        plt.xlabel('z')
        plt.grid(True)
        plt.ylim(0.002,0.005)
        for i in range(len(models)):
            da_distrib = da_list[i]
            z_array = np.ones(len(da_distrib))*1089
            plt.scatter(z_array, da_distrib, s=40, facecolors='none',
                        edgecolors=c[i], label=models[i])
        plt.legend()
        plt.show()

        plt.figure()
        plt.title(f'Histogram of angular diameter distances at z = {zpicks[-1]},'+
        f'\n $\sigma$ on data = {noise}, number of SN Ia used = {npoints}')
        plt.xlabel('$(H_0 /c) * D_A$')
        for i in range(len(models)):
            da_distrib = da_list[i]
            plt.hist(da_distrib, bins = 50, normed=True, color=c[i], histtype='step',
                     stacked=True, fill=False, label=models[i])
        plt.legend()
        plt.show()

        plt.figure()
        plt.title(f'$\mu$ and $\sigma$ of the angular diameter distance distribution'+
        f'\n $\sigma$ on data = {noise}, number of SN Ia used = {npoints}')
        plt.xlabel('$(H_0 /c) * D_A$')
        plt.ylabel(r'$z = 1089$')
        for i in range(len(models)):
            da_distrib = np.asarray(da_list[i])
            da_mean = np.mean(da_distrib)
            da_sd = np.std(da_distrib)
            if da_sd < da_mean/100:
                print(f'da_sd = {da_sd}, model = {test_key}, mean/100 = {da_mean/100}')
            plt.errorbar(da_mean, yaxis_tick[i], xerr=da_sd, fmt='-o', color=c[i],
                         ecolor=ec[i], elinewidth=3, capsize=0, label=models[i])
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