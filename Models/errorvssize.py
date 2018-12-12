#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 18:31:51 2018

@author: BallBlueMeercat
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
            ,'stepfall'
            ,'waterfall'
#            ,'exotic'           # NaN field at recombination
#            ,'late_intxde'
#            ,'heaviside_late_int'
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
    sigma = 0.001
    sigma_max = 0.2
    sigma_step = 0.7
    npoints_min = 104800
    npoints_max = 104801
    npoints_step = 1000
    if key:
        run = 0
        if key == 'waterfall':
            names = ['Mcorr',
                     'm_ombar', 'r_ombar', 'a_ombar', 'b_ombar', 'c_ombar',
                     'v_in', 'w_in', 'x_in', 'y_in', 'z_in']
            values = np.array([-19.3,
                               0.3, 0.025, 0.1, 0.1, 0.1,
                               0.0, 0.0, 0.0, 0.0, 0.0])
        elif key == 'stepfall':
            names = ['Mcorr', 'm_ombar', 'r_ombar', 'a_ombar',
                     'v_in', 'w_in', 'x_in']
            values = np.array([-19.3, 0.3, 0.025, 0.1, 0.0, 0.0, 0.0])
        elif key == 'exotic':
            names = ['Mcorr', 'm_ombar', 'r_ombar', 'gamma', 'zeta']
            values = np.array([-19.3, 0.3, 0.025, 0.0, 0.0])
        elif key == 'LCDM':
            names = ['Mcorr', 'm_ombar']
            values = np.array([-19.3, 0.3])
        else:
            names = ['Mcorr', 'm_ombar','gamma']
            values = np.array([-19.3, 0.3, 0.0])

        # Making sure number of parameters matches number of names given:
        assert len(names) == len(values), "len(names) != len(values)"

        # Folder for saving output.
        directory = f'{int(time.time())}_{key}'
        # Relative path of output folder.
        save_path = os.path.join('results_error_vs_data',directory)
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        while sigma < sigma_max:
            npoints = npoints_min
            while npoints < npoints_max:
                run += 1
                data_path = f'data/{npoints}_{all_zpicks[-1]}_sigma_{sigma}.p'
                my_file = Path(data_path)
                if my_file.is_file():
                    with open(data_path,'rb') as rfp: zpicks, nmag = pickle.load(rfp)
                else:
                    # Generating redshifts.
#                    n = len(zpicks)//npoints
#                    print(f'n = {n}')
#                    zpicks = all_zpicks[0::n]
#                    zpicks[1] = 0.0001
#                    zpicks[-1] = 3
#                    npoints = len(zpicks)
#                    print(f'npoints = {npoints}')
                    zpicks = all_zpicks
                    data_dic = {'zpicks':zpicks}
                    # Generating LCDM mag and da.
                    mag, da = datasim.magn(['Mcorr', 'matter'], np.array([-19.3, 0.3]), data_dic, 'LCDM', plot_key=False)
                    # Adding noise to LCDM mag.
                    nmag = datasim.gnoise(mag, mu, sigma)

                    plt.figure()
                    plt.title('Artificial data')
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
                npoints += npoints_step

            sigma += sigma_step
            sigma = round(sigma, 2)
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