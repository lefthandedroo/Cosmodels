#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 07:57:27 2019

@author: BallBlueMeercat

open propert, sampler from errorvssize runs and compare properties of resulting parmater distributions
"""

import matplotlib.pyplot as plt
import time
from pathlib import Path
import tools
import os.path
import pickle

import matplotlib as mpl
#mpl.style.use('default') # has to be switched on to set figure size
mpl.style.use('fivethirtyeight')
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['grid.color'] = 'white'

# Script timer.
timet0 = time.time()

#Model, errors on data and dataset sizes to iterate through:
test_keys = [None
#            ,'rainbow'
#            ,'kanangra'
#            ,'waterfall'
#            ,'stepfall'
#            ,'exotic'
#            ,'late_intxde'
#            ,'heaviside_late_int'
#            ,'heaviside_sudden'
#            ,'late_int'
#            ,'expgamma'
#            ,'txgamma'         # doesn't converge
#            ,'zxgamma'
#            ,'gamma_over_z'    # doesn't converge
#            ,'zxxgamma'        # gamma forced positive in firstderivs
#            ,'gammaxxz'        # gamma forced positive in firstderivs
#            ,'rdecay_m'
#            ,'rdecay_de'
#            ,'rdecay_mxde'
#            ,'rdecay'
#            ,'interacting'
            ,'LCDM'
#            ,'rLCDM'
            ]

sd_list = []
mean_list = []
vc_list = []

sigma_list = []
npoints_list = []
sampler_list = []

sigma_options = None, 0.01#, 0.07, 0.14, 0.2  #0.14, 0.2 #0.0001, 0.005, 0.007
npoints_options = 1048, 10480, 104800, 1048000

#N = len(sigma_options) * (len(npoints_options)) # Can't have None in options
N = len(npoints_options) # Ca't have None in options


for test_key in test_keys:
    if test_key:
        # parameter names and values used as input
        names, values = tools.names_values(test_key)
        # Relative path to folder for saving output.
        save_path = os.path.join('results_error_vs_data_plots',test_key)
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        for sigma in sigma_options:
            if sigma:
                pass
            else:
                continue
            for npoints in npoints_options:
                if npoints:
                    pass
                else:
                    continue
                #retrieving results from errovsdatasize runs
                file_path = f'results_error_vs_data/{test_key}/sigma{sigma}_npoints{npoints}.p'
                my_file = Path(file_path)
                if my_file.is_file():
                    # results of emcee runs
                    with open(file_path,'rb') as rfp: propert, sampler = pickle.load(rfp)
                    file_open_indicator = 1
#                    print(f'{file_path} opened successfully')
                else:
                    print(f'{file_path} failed to open')
                    file_open_indicator = 0
                # assert ensures only sampler from intended file_path is used
                assert file_open_indicator > 0, "file_open_indicator failed"
                for key in propert:
                    if 'sd' in key:
                        sd = propert.get(key,0)
                        sd_list.append([key, sd])
                    elif 'mean' in key:
                        mean = propert.get(key,0)
                        mean_list.append([key, mean])
                        if mean != 0:
                            vc = sd/mean * 100
                            vc_list.append([key[0]+'_vc', vc])

                sigma_list.append(sigma)
                npoints_list.append(npoints)
                sampler_list.append(sampler)

        for j in range(len(values)):
            sd = []
            mean = []
            vc = []
            for i in range(N):
                index = i*len(values)+j

                sd_name = sd_list[index][0]
                sd_initial = sd_name[0]
                sd.append(sd_list[index][1])

                mean_name = mean_list[index][0]
                mean_initial = mean_name[0]
                mean.append(mean_list[index][1])

                vc_name = vc_list[index][0]
                vc_initial = vc_name[0]
                vc.append(vc_list[index][1])

                i+=1

            fig, ax = plt.subplots()
            ax.scatter(npoints_list, sd, c="C{}".format(0))

            # Plotting sd vs dataset size.
            for i, txt in enumerate(sigma_list):
                txt = str(txt) # noise on data
                ax.annotate(txt, (npoints_list[i], sd[i]))

            plt.xlabel('Dataset size')
            plt.ylabel('$\sigma$')
            plt.title(sd_name+' vs dataset size'+
                      '\n s.d. of noise labeled, model '+test_key)
            stamp = str(int(time.time()))
            filename = str(stamp)+'_sd_of_'+sd_initial+'_.png'
            filename = os.path.join(save_path, filename)
            plt.savefig(filename)

            # Plotting mean vs dataset size.
            fig, ax = plt.subplots()
            ax.scatter(npoints_list, mean, c="C{}".format(1))
            for i, txt in enumerate(sigma_list):
                txt = str(txt) # noise on data
                ax.annotate(txt, (npoints_list[i], mean[i]))

            plt.xlabel('Dataset size')
            plt.ylabel('$\mu$')
            plt.title(mean_name+' vs dataset size'+
                      '\n s.d. of noise labeled, model '+test_key)
            stamp = str(int(time.time()))
            filename = str(stamp)+'_mean_of_'+mean_initial+'_.png'
            filename = os.path.join(save_path, filename)
            plt.savefig(filename)

            # Plotting variance coefficient vs dataset size.
            if len(vc) == N:
                fig, ax = plt.subplots()
                ax.scatter(npoints_list, vc, c="C{}".format(2))
                for i, txt in enumerate(sigma_list):
                    txt = str(txt)
#                    ax.annotate(txt, (npoints_list[i], vc[i]))

                plt.xlabel('Dataset size')
                plt.ylabel(r'$\sigma/ \mu \times 100$')
                plt.title(vc_name+' vs dataset size'+
                          '\n s.d. of noise labeled, model '+test_key)
                stamp = str(int(time.time()))
                filename = str(stamp)+'_vc_of_'+vc_initial+'_.png'
                filename = os.path.join(save_path, filename)
                plt.savefig(filename)

                j+=1

        plt.show()

        # Saving results to directory.
        import results
        results.save(save_path, 'vc_list', vc_list)
        results.save(save_path, 'sd_list', sd_list)
        results.save(save_path, 'mean_list', mean_list)

        results.save(save_path, 'sigma_list', sigma_list)
        results.save(save_path, 'npoints_list', npoints_list)
        results.save(save_path, 'sampler_list', sampler_list)

#        import plots
#        plots.precise_runs(key, names, values, 9, 1.2)

# Time taken by script.
timet1=time.time()
tools.timer('errorvsdatasize_distribs', timet0, timet1)