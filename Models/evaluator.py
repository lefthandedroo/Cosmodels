#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 21:45:40 2018

@author: BallBlueMeercat
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import time
import os.path
import results
import tools
import datasim
import stats
import plots


# Number of emcee steps.
nsteps = 100000

# Statistical parameteres of noise: mean, standard deviation.
mu, sigma = 0.0, 0.2 # sigma != 0

# Pantheon data:
dataname = './data/lcparam_full_long.txt'
pantheon = pd.read_csv(dataname, sep=" ")

# Reading each txt file column of interest as numpy.ndarray
mag = pantheon.mb.values
x1 = pantheon.x1.values
colour = pantheon.color.values
zpicks = pantheon.zhel.values

# Stacking arrays together and sorting by accending redshift.
data = np.stack((mag,x1,colour,zpicks), axis=0)
data.sort(axis=-1)

mag = data[0]
#x1 = data[1]
#colour = data[2]
zpicks = data[3]

# Pantheon data.
data_dic = {'mag':mag, 'zpicks':zpicks}

## Generating redshifts and LCDM mag and da.
#zpicks = np.random.uniform(low=0.0, high=1101, size=(1000,))
#zpicks = np.sort(zpicks, axis=None)
#data_dic = {'mag':None, 'zpicks':zpicks}
#names = ['Mcorr', 'matter']
#values = np.array([-19.3, 0.3])
#mag = datasim.noisy_mag(mu, sigma, names, values, data_dic, 'LCDM')
#data_dic = {'mag':mag, 'zpicks':zpicks}

## Plot param evolutions for multiple models on the same axis.
#p1 = ['Mcorr', 'm_ombar'], np.array([-19.3, 0.0])
#p2 = ['Mcorr', 'm_ombar', 'r_ombar', 'a_ombar','v_in', 'w_in', 'x_in'], np.array([-19.3, 0.3, 0.025, 0.1, 0.0, 0.0, 0.0])
#datasim.model_comparison([p1, p2], zpicks, ['LCDM', 'stepfall'], plot_key=True)


firstderivs_functions = [None
            ,'stepfall'
            ,'waterfall'
            ,'exotic'
            ,'late_intxde'
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

def modelcheck():

    for test_key in firstderivs_functions:
        if test_key:
            print('---',test_key)
            if test_key == 'waterfall':
                names = ['Mcorr',
                         'm_ombar', 'r_ombar', 'a_ombar', 'b_ombar', 'c_ombar',
                         'v_in', 'w_in', 'x_in', 'y_in', 'z_in']
                values = np.array([-19.3,
                                   0.3, 0.025, 0.1, 0.1, 0.1,
                                   0.0, 0.0, 0.0, 0.0, 0.0])
            elif test_key == 'stepfall':
                names = ['Mcorr', 'm_ombar', 'r_ombar', 'a_ombar',
                         'v_in', 'w_in', 'x_in']
                values = np.array([-19.3, 0.3, 0.025, 0.1, 0.0, 0.0, 0.0])
            elif test_key == 'exotic':
                names = ['Mcorr', 'm_ombar', 'r_ombar', 'gamma', 'zeta']
                values = np.array([-19.3, 0.3, 0.025, 0.0, 0.0])
            elif test_key == 'LCDM':
                names = ['Mcorr', 'm_ombar']
                values = np.array([-19.3, 0.3])
            else:
                names = ['Mcorr', 'm_ombar','gamma']
                values = np.array([-19.3, 0.3, 0.0])

            # Making sure number of parameters matches number of names given:
            if len(names) != len(values):
                raise ValueError('len(names) != len(values)')
            datasim.magn(names, values, data_dic, test_key, plot_key=True)
    return

#modelcheck()

def emcee():
    print('@@@@@@@ Mcor_emcee @@@@@@@')

    for test_key in firstderivs_functions:

        if test_key:
            print('---',test_key)

            if test_key == 'waterfall':
                names = ['Mcorr',
                         'm_ombar', 'r_ombar', 'a_ombar', 'b_ombar', 'c_ombar',
                         'v_in', 'w_in', 'x_in', 'y_in', 'z_in']
                values = np.array([-19.3,
                                   0.3, 0.025, 0.1, 0.1, 0.1,
                                   0.0, 0.0, 0.0, 0.0, 0.0])
            elif test_key == 'stepfall':
                names = ['Mcorr', 'm_ombar', 'r_ombar', 'a_ombar',
                         'v_in', 'w_in', 'x_in']
                values = np.array([-19.3, 0.3, 0.025, 0.1, 0.0, 0.0, 0.0])
            elif test_key == 'exotic':
                names = ['Mcorr', 'm_ombar', 'r_ombar', 'gamma', 'zeta']
                values = np.array([-19.3, 0.3, 0.025, 0.0, 0.0])
            elif test_key == 'LCDM':
                names = ['Mcorr', 'm_ombar']
                values = np.array([-19.3, 0.3])

            else:
                names = ['Mcorr', 'm_ombar','gamma']
                values = np.array([-19.3, 0.3, 0.0])

            # Making sure number of parameters matches number of names given:
            if len(names) != len(values):
                raise ValueError('len of parameters names and values are not matching')

            # Creating a folder for saving output.
            save_path = './results_emcee/'+str(int(time.time()))+'_'+test_key
            if not os.path.exists(save_path):
                os.makedirs(save_path)

            # Script timer.
            timet0 = time.time()

            # emcee parameter search.
            propert, sampler = stats.stats(names, values, data_dic, sigma,
                                           nsteps, save_path, test_key, plot=1)
            # Time taken by script.
            timet1=time.time()
            tools.timer('script', timet0, timet1)

            # Saving sampler to directory.
            results.save(save_path, 'sampler', sampler)

            print('Data:',dataname)

    return

emcee()

def errorvsdatasize():

    test_key = 'exotic'
    params_dic = [{'matter':0.3},{'Mcorr':-19.3}, {'gamma':0.0}, {'zeta':0.0}]

    # Script timer.
    timet0 = time.time()

    sigma = 0.3
    sigma_max = 1
    sigma_step = 0.2
    npoints_min = 2750
    npoints_max = 4000
    npoints_step = 600

    # How many iterations have I signed up for?
    N = tools.runcount(sigma, sigma_max, sigma_step,
                   npoints_min, npoints_max, npoints_step)

    if input('Happy with '+str(N)+' iterations? (enter=yes) '):
        return

    # Folder for saving output.
    directory = str(int(time.time()))+'_'+test_key
    # Relative path of output folder.
    save_path = './results_error_vs_data/'+test_key+'/'+directory
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    run = 0

    sd_list = []
    mean_list = []
    vc_list = []

    sigma_list = []
    npoints_list = []
    sampler_list = []

    while sigma < sigma_max:

        npoints = npoints_min

        while npoints < npoints_max:
            print('_____________________ run number',run)

            zpicks = datasim.redshift_picks(0.005, 2.0, npoints)
            data_dic = {'zpicks':zpicks}

            # Magnitudes corresponding to simulated redshifts:
            mag = datasim.noisy_mag(mu, sigma, params_dic, data_dic, 'LCDM')

            data_dic['mag']=mag

            propert, sampler = stats.stats(params_dic, data_dic,
                                           sigma, nsteps, save_path,
                                           test_key, plot=False)

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

            npoints += npoints_step
            run += 1

        sigma += sigma_step
        sigma = round(sigma, 2)

    for j in range(len(params_dic)):
        sd = []
        mean = []
        vc = []
        for i in range(N):
            index = i*len(params_dic)+j

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
        ax.scatter(npoints_list, sd, c='r')

        # Plotting SD vs dataset size.
        for i, txt in enumerate(sigma_list):
            txt = 'sd = '+ str(txt)
            ax.annotate(txt, (npoints_list[i], sd[i]))

        plt.xlabel('Dataset size')
        plt.ylabel('s.d. of a marginalised distribution')
        plt.title(sd_name+' vs dataset size'+
                  '\n s.d. of noise labeled, model '+test_key)
        stamp = str(int(time.time()))
        filename = str(stamp)+'_sd_of_'+sd_initial+'_.png'
        filename = os.path.join(save_path, filename)
        plt.savefig(filename)

        # Plotting mean vs dataset size.
        fig, ax = plt.subplots()
        ax.scatter(npoints_list, mean, c='c')
        for i, txt in enumerate(sigma_list):
            txt = 'sd = '+ str(txt)
            ax.annotate(txt, (npoints_list[i], mean[i]))

        plt.xlabel('Dataset size')
        plt.ylabel('Mean of a marginalised distribution')
        plt.title(mean_name+' vs dataset size'+
                  '\n s.d. of noise labeled, model '+test_key)
        stamp = str(int(time.time()))
        filename = str(stamp)+'_mean_of_'+mean_initial+'_.png'
        filename = os.path.join(save_path, filename)
        plt.savefig(filename)

        # Plotting variance coefficient vs dataset size.
        if len(vc) == N:
            fig, ax = plt.subplots()
            ax.scatter(npoints_list, vc, c='g')
            for i, txt in enumerate(sigma_list):
                txt = 'sd = '+ str(txt)
                ax.annotate(txt, (npoints_list[i], vc[i]))

            plt.xlabel('Dataset size')
            plt.ylabel('s.d. /mean x100 of a marginalised distribution')
            plt.title(vc_name+' vs dataset size'+
                      '\n s.d. of noise labeled, model '+test_key)
            stamp = str(int(time.time()))
            filename = str(stamp)+'_vc_of_'+vc_initial+'_.png'
            filename = os.path.join(save_path, filename)
            plt.savefig(filename)

            j+=1

    plt.show()

#     Saving results to directory.
    results.save(save_path, 'vc_list', vc_list)
    results.save(save_path, 'sd_list', sd_list)
    results.save(save_path, 'mean_list', mean_list)

    results.save(save_path, 'sigma_list', sigma_list)
    results.save(save_path, 'npoints_list', npoints_list)
    results.save(save_path, 'sampler_list', sampler_list)

    print('directory:',directory)

    # Time taken by script.
    timet1=time.time()
    tools.timer('errorvsdatasize', timet0, timet1)

    plots.precise_runs(test_key, params_dic, 9, 1.2)

    return vc_list, sd_list, mean_list, sigma_list, npoints_list, sampler_list

#vc, sd, mean, sigma, npoints, sampler = errorvsdatasize()

