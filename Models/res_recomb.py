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

import matplotlib as mpl
#mpl.style.use('default') # has to be switched on to set figure size
mpl.style.use('fivethirtyeight')
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['grid.color'] = 'white'

zpicks = np.array([1089])
data_dic = {'zpicks':zpicks}

timed = False

if timed:
    import cProfile, pstats, io
    pr = cProfile.Profile()
    pr.enable()

models = 'exotic', 'LCDM'
noise_options = 0.01, None      #0.001, 0.01, 0.07, 0.14, 0.2
npoints_options = 1048, 10480, 1048000  #1048, 10480, 104800, 1048000
yaxis_tick = 1, 2, 3 #fake y values for visual separation of da's. da's plotted are all at z=1089.
msize = 70, 20  # marker sizes to differentiate between scatter plots
n_bin = 100 # histogram bin number

for npoints in npoints_options:
    if not npoints: # skipping the None to allow comparing a list of 1
        continue
    for noise in noise_options:
        if not noise:
            continue
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
                # get da for every skip-th set of params, keep skip < n_bin
                # skip=50 and n_bin=50 works
                skip = 10
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

        plt.figure()
        plt.title(f'$\sigma$ on data = {noise}, {npoints} SN Ia used')
        plt.ylabel(r'$z = 1089$')
        plt.xlabel('$(H_0 /c) * D_A$')
        plt.grid(True)
        plt.xlim(0.00285,0.00295)
        for i in range(len(models)):
            da_distrib = da_list[i]
            z_array = np.ones(len(da_distrib))*1089
            plt.scatter(da_distrib, z_array, s=msize[i], facecolors='none', edgecolors="C{}".format(i), label=models[i])
        plt.locator_params(axis='x', nbins=4)
        plt.yticks([])
        plt.legend()
        plt.show()

#        # Non-normalised histogram
#        plt.figure()
#        plt.title(f'$D_A$ at z = {zpicks[-1]},'
#                +f'\n $\sigma$ on data = {noise}, {npoints} SN Ia used')
#        plt.xlabel('$(H_0 /c) * D_A$')
#        for i in range(len(models)):
#            da_distrib = da_list[i]
#            plt.hist(da_distrib, normed=False, histtype='step',
#                     stacked=True, label=models[i])
#        plt.locator_params(axis='x', nbins=5)
#        plt.legend()
#        plt.show()

        cutoff = 0.1
        plt.figure()
        smallest_da = 1
        largest_da = 0
        plt.title(f'scatter of $D_A$ histogram $y > {cutoff}$'
                  +'\n from all guessed parameter sets')
        for i in range(len(models)):
            face_color = 'none', "C{}".format(i)
            da_distrib = da_list[i]
            y, x = np.histogram(da_distrib, bins=n_bin)
#            print('x = ',x)
            y_norm = y/max(y)
            x = x[1:]
#            print('scatter x[1:] = ',x)
#            print('scatter y = ',y)
#            print('scatter y_norm = ',y_norm)
            indices = np.where(y_norm > cutoff)
            y_norm = y_norm[indices]
            x = x[indices]
            if min(x) < smallest_da:
                smallest_da = min(x)
            if max(x) > largest_da:
                largest_da = max(x)
            plt.scatter(x, y_norm, s=msize[i], facecolors=face_color[i], edgecolors="C{}".format(i), label=f'{models[i]}')
        plt.xlim((smallest_da-0.000001),(largest_da+0.000001))
        plt.locator_params(axis='x', nbins=5)
        plt.legend()
        plt.show()


        plt.figure()
        smallest_da = 1
        largest_da = 0
        plt.title(f'scatter of $D_A$ histogram $y > {cutoff}$'
                  +'\n from all guessed parameter sets')
        for i in range(len(models)):
            face_color = 'none', "C{}".format(i)
            da_distrib = da_list[i]
            y, x = np.histogram(da_distrib, bins=n_bin)
#            print('x = ',x)
            y_norm = y/max(y)
            x = x[:-1]
#            print('scatter x[:-1] = ',x)
#            print('scatter y = ',y)
#            print('scatter y_norm = ',y_norm)
            indices = np.where(y_norm > cutoff)
            y_norm = y_norm[indices]
            x = x[indices]
            if min(x) < smallest_da:
                smallest_da = min(x)
            if max(x) > largest_da:
                largest_da = max(x)
            plt.scatter(x, y_norm, s=msize[i], facecolors=face_color[i], edgecolors="C{}".format(i), label=f'{models[i]}')
        plt.xlim((smallest_da-0.000001),(largest_da+0.000001))
        plt.locator_params(axis='x', nbins=5)
        plt.legend()
        plt.show()

        plt.figure()
        plt.title('scatter of complete $D_A$ histogram'
                  +'\n from all guessed parameter sets')
        for i in range(len(models)):
            face_color = 'none', "C{}".format(i)
            da_distrib = da_list[i]
            y, x = np.histogram(da_distrib, bins=n_bin)
#            print('scatter x = ',x)
            y_norm = y/max(y)
            x = x[:-1]
#            print('scatter x[:-1] = ',x)
#            print('scatter y = ',y)
#            print('scatter y_norm = ',y_norm)
            plt.scatter(x, y_norm, s=msize[i], facecolors=face_color[i], edgecolors="C{}".format(i), label=f'{models[i]}')
        plt.xlim(0.00284,0.00294)
        plt.locator_params(axis='x', nbins=5)
        plt.legend()
        plt.show()





#        for i in range(len(models)):
#            plt.figure()
#            da_distrib = da_list[i]
#            y, x = np.histogram(da_distrib)
#            y_norm = y/max(y)
#            x = x[1:]
#            plt.title("$D_A$'s using all guessed parameters")
#            plt.plot(x, y_norm, label=f'{models[i]}')
#            plt.xlim(0.00285,0.003)
#            plt.locator_params(axis='x', nbins=5)
#            plt.legend()
#        plt.show()

#        plt.figure()
#        plt.title(f'$\mu$ and $\sigma$ of the $D_A$ distribution'
#                  +f'\n $\sigma$ on data = {noise}, {npoints} SN Ia used')
#        plt.xlabel('$(H_0 /c) * D_A$')
#        plt.ylabel(r'$z = 1089$')
#        for i in range(len(models)):
#            da_distrib = np.asarray(da_list[i])
#            da_mean = np.mean(da_distrib)
#            da_sd = np.std(da_distrib)
#            if da_sd < da_mean/100:
#                print(f'{test_key} da_sd (={da_sd}) < mean/100 (={da_mean/100})')
##            plt.errorbar(da_mean, yaxis_tick[i], xerr=da_sd, fmt='-o',
##        color=c[i], ecolor=ec[i], elinewidth=3, capsize=0, label=models[i])
#            plt.errorbar(da_mean, yaxis_tick[i], xerr=da_sd, fmt='o', label=models[i])
#        plt.locator_params(axis='x', nbins=5)
#        plt.ylim(-4,8)
#        plt.yticks([])
#        plt.legend()
#        plt.show()


if timed:
    pr.disable()
    s = io.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print (s.getvalue())