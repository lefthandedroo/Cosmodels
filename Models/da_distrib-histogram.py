#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 11:45:04 2019

@author: BallBlueMeercat

This script runs through models, noise on data and dataset size options
to collect the distribution of angular diameter distance assosiated with each,
then plots these distributions as scatter and as a normalised histogram.
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.style.use('default') # has to be switched on to set figure size
mpl.style.use('fivethirtyeight')
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['grid.color'] = 'white'

timed = False

if timed:
    import cProfile, pstats, io
    pr = cProfile.Profile()
    pr.enable()

models = 'LCDM', 'exotic', 'stepfall'
noise_options = 0.001, None      #0.001, 0.01, 0.07, 0.14, 0.2
npoints_options = None, 104800  #1048, 10480, 104800, 1048000
#yaxis_tick = fake y values for visual separation of da's (which are all at z=1089)
yaxis_tick = 1, 2, 3
msize = 70, 20, 10  # marker sizes to differentiate between scatter plots
n_bin = 30, 10000 # histogram bin number
n_bin = 30, 10000, 10000
f_dic = {}  # all non-normalised frequencies of D_A, as well as non-normalised frequencies corresponding to indecies where normalised y > 0.1.

for npoints in npoints_options:
    if not npoints: # skipping None to allow iteratng over a list of 1
        continue
    for noise in noise_options:
        if not noise:
            continue
        da_list = []
        for test_key in models:
            if test_key:
                filename = f'da_distribs/da_distrib_{test_key}_{noise}_{npoints}.p'
                try:
                    with open(filename,'rb') as rfp: da_distrib = pickle.load(rfp)
                    file_open_indicator = 1
                except:
                    print("didnt't open",filename)
                    file_open_indicator = 0
                # assert ensures only da_distrib from intended file_path is used
                assert file_open_indicator > 0, f"{filename} didn't have a da_distrib"
                da_list.append(da_distrib)

        # Scatter plot of all D_A
        smallest_da = 1
        largest_da = 0
        plt.figure()
#        plt.title(f'$\sigma$ on data = {noise}, {npoints} SN Ia used')
        plt.ylabel(r'$z = 1089$')
        plt.xlabel('$D_A (H_0 /c)$')
        plt.grid(True)
        for i in range(len(models)):
            da_distrib = da_list[i]
            if min(da_distrib) < smallest_da:
                smallest_da = min(da_distrib)
            if max(da_distrib) > largest_da:
                largest_da = max(da_distrib)
            z_array = np.ones(len(da_distrib))*1089
            plt.scatter(da_distrib, z_array, s=msize[i], facecolors='none', edgecolors="C{}".format(i), label=models[i])
        plt.xlim((smallest_da-0.0000015),(largest_da+0.0000015))
        plt.locator_params(axis='x', nbins=4)
        plt.yticks([])
        plt.legend()
        plt.show()

        # Scatter of a normalised x[1:] D_A histogram, only values >0.1
        cutoff = 0.1
        plt.figure()
        smallest_da = 1
        largest_da = 0
#        plt.title(f'scatter of $D_A$ histogram $y > {cutoff}$'
#                  +'\n from all guessed parameter sets')
        plt.xlabel('$D_A (H_0 /c)$')
        plt.ylabel('$f$')
        for i in range(len(models)):
            face_color = 'none', "C{}".format(i), "C{}".format(i)
            da_distrib = da_list[i]
            y, x = np.histogram(da_distrib, bins=n_bin[i])
            y_norm = y/max(y)
            x = x[1:]
            indices = np.where(y_norm > cutoff)
            f_dic[models[i]+' > 0.1'] = y[indices]
            f_dic[models[i]+' all'] = y
            y_norm = y_norm[indices]
            x = x[indices]
            if min(x) < smallest_da:
                smallest_da = min(x)
            if max(x) > largest_da:
                largest_da = max(x)
            print(f'{models[i]} x= {x}')
            plt.scatter(x, y_norm, s=msize[i], facecolors=face_color[i], edgecolors="C{}".format(i), label=f'{models[i]}')
        plt.xlim((smallest_da),(largest_da))
        plt.locator_params(axis='x', nbins=5)
        plt.legend()
        plt.show()

#        # Scatter of a normalised x[:-1] D_A histogram, only values >0.1
#        plt.figure()
#        smallest_da = 1
#        largest_da = 0
#        plt.title(f'scatter of x[:-1] $D_A$ histogram $y > {cutoff}$'
#                  +'\n from all guessed parameter sets')
#        for i in range(len(models)):
#            face_color = 'none', "C{}".format(i)
#            da_distrib = da_list[i]
#            y, x = np.histogram(da_distrib, bins=n_bin[i])
#            y_norm = y/max(y)
#            x = x[:-1]
#            indices = np.where(y_norm > cutoff)
#            y_norm = y_norm[indices]
#            x = x[indices]
#            if min(x) < smallest_da:
#                smallest_da = min(x)
#            if max(x) > largest_da:
#                largest_da = max(x)
#            plt.scatter(x, y_norm, s=msize[i], facecolors=face_color[i], edgecolors="C{}".format(i), label=f'{models[i]}')
#        plt.xlim((smallest_da-0.000001),(largest_da+0.000001))
#        plt.locator_params(axis='x', nbins=5)
#        plt.legend()
#        plt.show()

        # Scatter of a normalised x[:-1] D_A histogram, all values.
        plt.figure()
        smallest_da = 1
        largest_da = 0
#        plt.title('scatter of complete $D_A$ histogram'
#                  +'\n from all guessed parameter sets')
        plt.xlabel('$D_A (H_0 /c)$')
        plt.ylabel('$f$')
        for i in range(len(models)):
            face_color = 'none', "C{}".format(i), "C{}".format(i)
            da_distrib = da_list[i]
            y, x = np.histogram(da_distrib, bins=n_bin[i])
            y_norm = y/max(y)
            x = x[:-1]
            if min(x) < smallest_da:
                smallest_da = min(x)
            if max(x) > largest_da:
                largest_da = max(x)
            plt.scatter(x, y_norm, s=msize[i], facecolors=face_color[i], edgecolors="C{}".format(i), label=f'{models[i]}')
        plt.xlim((smallest_da),(largest_da))
        plt.locator_params(axis='x', nbins=5)
        plt.legend()
        plt.show()

if timed:
    pr.disable()
    s = io.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print (s.getvalue())