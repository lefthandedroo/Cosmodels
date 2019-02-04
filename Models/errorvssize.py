#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 18:31:51 2018

@author: BallBlueMeercat

Runs through given sigma (error on data) and npoints (dataset size), use pre-made artificial mag or create new mag using pre-made zpicks. Saves sampler.p for each run.
"""
import matplotlib.pyplot as plt
import numpy as np
import time
import datasim
from pathlib import Path
import tools
import os.path
import stats
import pickle

# Script timer.
timet0 = time.time()

nsteps = 10000   # Number of emcee steps.
mu = 0.0        # Mean of noise added to LCDM to simulate data.

#Model, errors on data and dataset sizes to iterate through:
test_keys = [None
#            ,'rainbow'
#            ,'kanangra'
#            ,'waterfall'
#            ,'stepfall'
            ,'exotic'
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

# Extracting pre-made redshifts z=0 to z=3.
try:
    with open('data/zpicks_1048_3.p','rb') as rfp: all_zpicks = pickle.load(rfp)
except:
    print("zpicks_1048_3.p didnt't open")

#sd_list = []
#mean_list = []
#vc_list = []
#
#sigma_list = []
#npoints_list = []
#sampler_list = []

for key in test_keys:
    sigma_options = 0.01, 0.07 #0.14, 0.2 #0.0001, 0.005
    npoints_options = 1048000, 104800#1048, 10480, 104800 #1048000, 10480000
    if key:
        run = 0
        names, values = tools.names_values(key)
        # Folder for saving output.
        directory = f'{int(time.time())}_{key}'
        # Relative path of output folder.
        save_path = os.path.join('results_error_vs_data',directory)
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        for sigma in sigma_options:
            for npoints in npoints_options:
                run += 1
                data_path = f'data/{npoints}_{all_zpicks[-1]}_sigma_{sigma}.p'
                my_file = Path(data_path)
                if my_file.is_file():
                    with open(data_path,'rb') as rfp: zpicks, nmag = pickle.load(rfp)
                else:
                    print(f'failed to get zpicks, nmag from {data_path}')
                    # Generating redshifts.
                    zpicks = np.random.uniform(low=0.0001, high=3, size=(npoints,))
                    zpicks = np.sort(zpicks, axis=None)
                    if zpicks[-1] != 3:
                        zpicks[-1] = 3
                    data_dic = {'zpicks':zpicks}
                    # Generating LCDM mag and da.
                    mag, da = datasim.magn(['Mcorr', 'matter'], np.array([-19.3, 0.3]), data_dic, 'LCDM', plot_key=False)
                    # Adding noise to LCDM mag.
                    nmag = datasim.gnoise(mag, mu, sigma)

                    plt.figure()
                    plt.title(f'Artificial data N={len(zpicks)}, $\sigma$={sigma}')
                    plt.scatter(zpicks, nmag)
                    plt.show()

                    data = zpicks, nmag
                    pickle.dump(data, open(data_path, 'wb'))

                data_dic = {'mag':nmag, 'zpicks':zpicks}

                print(f'--- {key} --------- run number {run}')
                propert, sampler = stats.stats(names, values, data_dic,
                                               sigma, nsteps, save_path,
                                               key, plot=False)
                output = propert, sampler
                output_path = os.path.join(save_path, f'sigma{sigma}_npoints{npoints}.p')
                pickle.dump(output, open(output_path, 'wb'))

#                import results
#                import plots
#                for key in propert:
#                    if 'sd' in key:
#                        sd = propert.get(key,0)
#                        sd_list.append([key, sd])
#                    elif 'mean' in key:
#                        mean = propert.get(key,0)
#                        mean_list.append([key, mean])
#                        if mean != 0:
#                            vc = sd/mean * 100
#                            vc_list.append([key[0]+'_vc', vc])
#
#                sigma_list.append(sigma)
#                npoints_list.append(npoints)
#                sampler_list.append(sampler)
#
#        for j in range(len(values)):
#            sd = []
#            mean = []
#            vc = []
#            for i in range(N):
#                index = i*len(values)+j
#
#                sd_name = sd_list[index][0]
#                sd_initial = sd_name[0]
#                sd.append(sd_list[index][1])
#
#                mean_name = mean_list[index][0]
#                mean_initial = mean_name[0]
#                mean.append(mean_list[index][1])
#
#                vc_name = vc_list[index][0]
#                vc_initial = vc_name[0]
#                vc.append(vc_list[index][1])
#
#                i+=1
#
#            fig, ax = plt.subplots()
#            ax.scatter(npoints_list, sd, c='r')
#
#            # Plotting SD vs dataset size.
#            for i, txt in enumerate(sigma_list):
#                txt = 'sd = '+ str(txt)
#                ax.annotate(txt, (npoints_list[i], sd[i]))
#
#            plt.xlabel('Dataset size')
#            plt.ylabel('s.d. of a marginalised distribution')
#            plt.title(sd_name+' vs dataset size'+
#                      '\n s.d. of noise labeled, model '+key)
#            stamp = str(int(time.time()))
#            filename = str(stamp)+'_sd_of_'+sd_initial+'_.png'
#            filename = os.path.join(save_path, filename)
#            plt.savefig(filename)
#
#            # Plotting mean vs dataset size.
#            fig, ax = plt.subplots()
#            ax.scatter(npoints_list, mean, c='c')
#            for i, txt in enumerate(sigma_list):
#                txt = 'sd = '+ str(txt)
#                ax.annotate(txt, (npoints_list[i], mean[i]))
#
#            plt.xlabel('Dataset size')
#            plt.ylabel('Mean of a marginalised distribution')
#            plt.title(mean_name+' vs dataset size'+
#                      '\n s.d. of noise labeled, model '+key)
#            stamp = str(int(time.time()))
#            filename = str(stamp)+'_mean_of_'+mean_initial+'_.png'
#            filename = os.path.join(save_path, filename)
#            plt.savefig(filename)
#
#            # Plotting variance coefficient vs dataset size.
#            if len(vc) == N:
#                fig, ax = plt.subplots()
#                ax.scatter(npoints_list, vc, c='g')
#                for i, txt in enumerate(sigma_list):
#                    txt = 'sd = '+ str(txt)
#                    ax.annotate(txt, (npoints_list[i], vc[i]))
#
#                plt.xlabel('Dataset size')
#                plt.ylabel('s.d. /mean x100 of a marginalised distribution')
#                plt.title(vc_name+' vs dataset size'+
#                          '\n s.d. of noise labeled, model '+key)
#                stamp = str(int(time.time()))
#                filename = str(stamp)+'_vc_of_'+vc_initial+'_.png'
#                filename = os.path.join(save_path, filename)
#                plt.savefig(filename)
#
#                j+=1
#
#        plt.show()
#
##     Saving results to directory.
#results.save(save_path, 'vc_list', vc_list)
#results.save(save_path, 'sd_list', sd_list)
#results.save(save_path, 'mean_list', mean_list)
#
#results.save(save_path, 'sigma_list', sigma_list)
#results.save(save_path, 'npoints_list', npoints_list)
#results.save(save_path, 'sampler_list', sampler_list)

print('directory:',directory)

# Time taken by script.
timet1=time.time()
tools.timer('errorvsdatasize', timet0, timet1)

#plots.precise_runs(key, names, values, 9, 1.2)