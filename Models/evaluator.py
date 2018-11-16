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
import pickle
import results
import tools
import datasim
import stats
import plots


# Number of emcee steps.
nsteps = 1000

# Statistical parameteres of noise: mean, standard deviation.
mu, sigma = 0.0, 0.07 # sigma != 0

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
zpicks = zpicks.tolist()

# Pantheon data.
data_dic = {'mag':mag, 'zpicks':zpicks}
#plt.figure()
#plt.title('Pantheon')
#plt.ylabel('Mag')
#plt.xlabel('redshift')
#plt.scatter(zpicks, mag)

# Plotting one model.
#zpicks[-6] = 5
#zpicks[-5] = 6
#zpicks[-4] = 7
#zpicks[-3] = 8
#zpicks[-2] = 9
#zpicks[-1] = 10

## Creating LCDM data.
#names = ['Mcorr', 'matter']
#values = np.array([-19.3, 0.3])
#mag = datasim.noisy_mag(mu, sigma, names, values, data_dic, 'LCDM')
#data_dic = {'mag':mag, 'zpicks':zpicks}

# Mag for LCDM.
names = ['Mcorr', 'matter']
values = np.array([-19.3, 0.3])
mag0 = datasim.magn(names, values, data_dic, 'LCDM', plot_key=False)

# Mag for waterfall in LCDM mode.
names = ['Mcorr','matter', 'radiation', 'a_ombar', 'b_ombar', 'c_ombar',
         'v_in', 'w_in', 'x_in', 'y_in', 'z_in']
values = np.array([-19.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
mag1 = datasim.magn(names, values, data_dic, 'waterfall', plot_key=False)

# Mag for waterfall wirg fluids but no interaction.
values = np.array([-19.3,0.3, 0.025, 0.1, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0])
mag2 = datasim.magn(names, values, data_dic, 'waterfall', plot_key=False)

# Sampler from chosen waterfall run.
with open('results_emcee/1541470181_waterfall/sampler.p','rb') as rfp:
    sampler = pickle.load(rfp)

# Mag from parameters with highest likelihood.
bi = np.argmax(sampler.flatlnprobability)
values = sampler.flatchain[bi,:]
mag3 = datasim.magn(names, values, data_dic, 'waterfall', plot_key=False)

# Parameters with x best likelihood.
flatlnprobability = sampler.flatlnprobability
flatchain_M = sampler.flatchain[:,0]
flatchain_m = sampler.flatchain[:,1]
flatchain_r = sampler.flatchain[:,2]
flatchain_a = sampler.flatchain[:,3]
flatchain_b = sampler.flatchain[:,4]
flatchain_c = sampler.flatchain[:,5]
flatchain_v = sampler.flatchain[:,6]
flatchain_w = sampler.flatchain[:,7]
flatchain_x = sampler.flatchain[:,8]
flatchain_y = sampler.flatchain[:,9]
flatchain_z = sampler.flatchain[:,10]
# Stacking them together and sorting by accending redshift.
flat_sorted = np.stack((flatchain_M, flatchain_m, flatchain_r,
                        flatchain_a, flatchain_b, flatchain_c,
                        flatchain_v, flatchain_w, flatchain_x,
                        flatchain_y, flatchain_z,flatlnprobability), axis=0)
flat_sorted.sort(axis=-1)

# Mag from parameters with second highest likelihood.
second_best=[]
for i in range(len(values)):
    second_best.append(sampler.flatchain[2,i])
values = np.asarray(second_best)
mag4 = datasim.magn(names, values, data_dic, 'waterfall', plot_key=False)

# Mag from parameters with lowest likelihood.
worst=[]
for i in range(len(values)):
    worst.append(sampler.flatchain[-2,i])
values = np.asarray(worst)
mag5 = datasim.magn(names, values, data_dic, 'waterfall', plot_key=False)

#plt.figure()
#plt.title('SN Ia magnitudes')
#plt.ylabel('Mag')
#plt.xlabel('redshift')
#plt.scatter(zpicks, mag, label='Pantheon set', marker=',', s=1)
#plt.plot(zpicks, mag0, label='LCDM')
#plt.plot(zpicks, mag1, label='waterfall in LCDM mode')
#plt.plot(zpicks, mag2, label='waterfall with fluids')
#plt.plot(zpicks, mag3, label='waterfall highest likelihood')
#plt.plot(zpicks, mag4, label='waterfall second highest likelihood')
#plt.plot(zpicks, mag5, label='waterfall lowest likelihood')
#plt.legend()

data_diff = mag - mag0
best_diff = mag0 - mag3
second_best_diff = mag0 - mag4
worst_diff = mag0 - mag5

#plt.figure()
#plt.title('Residuals')
#plt.scatter(zpicks, data_diff, label='data-LCDM', marker=',', s=1)
#plt.plot(zpicks, best_diff, label='LCDM - best emcee fit')
#plt.plot(zpicks, second_best_diff, label='LCDM - 2nd best emcee fit')
#plt.plot(zpicks, worst_diff, label='LCDM - worst emcee fit')
#plt.legend()



#names = ['Mcorr','matter', 'radiation','gamma', 'zeta']
#values = np.array([-19.3,0.3, 0.025, 0.0, 0.0])
#datasim.magn_plot(names, values, data_dic, 'exotic')
#
#names = ['Mcorr', 'matter']
#values = np.array([-19.3, 0.3])
#datasim.magn_plot(names, values, data_dic, 'LCDM', True)



# Compare param evolution for 3 models, plotting on the same axis.
#g2, g3, z3 = 0.0, 0.0, 0.0
#p1 = ['Mcorr', 'matter',], np.array([-19.3, 0.0])
#p2 = ['Mcorr', 'matter','gamma'], np.array([-19.3, 0.3, g2])
#p3 = ['Mcorr', 'matter', 'radiation', 'gamma', 'zeta'], np.array([-19.3, 0.3, 0.0, g3, z3])
#datasim.model_comparison([p1, p2, p3], zpicks, ['LCDM', 'late_int', 'exotic'],
#    ['no interaction','$\gamma$='+str(g2),'$\gamma$='+str(g3)+' $\zeta$='+str(z3)])

firstderivs_functions = [None
            ,'waterfall'
#            ,'exotic'
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

def modelcheck():

    for test_key in firstderivs_functions:
        if test_key:
            print('---',test_key)
            if test_key == 'waterfall':
                names = ['Mcorr',
                         'matter', 'radiation', 'a_ombar', 'b_ombar', 'c_ombar',
                         'v_in', 'w_in', 'x_in', 'y_in', 'z_in']
                values = np.array([-19.3,
                                   0.3, 0.025, 0.1, 0.1, 0.1,
                                   0.0, 0.0, 0.0, 0.0, 0.0])
            elif test_key == 'exotic':
                names = ['Mcorr', 'matter', 'radiation', 'gamma', 'zeta']
                values = np.array([-19.3, 0.3, 0.025, 0.0, 0.0])
            elif test_key == 'LCDM' or test_key == 'dlLCDM':
                names = ['Mcorr', 'matter']
                values = np.array([-19.3, 0.3])
            else:
                names = ['Mcorr', 'matter','gamma']
                values = np.array([-19.3, 0.3, 0.0])

            # Making sure number of parameters matches number of names given:
            if len(names) != len(values):
                raise ValueError('len(names) != len(values)')
            datasim.magn(names, values, data_dic, test_key, plot_key=True)
    return 

modelcheck()

def emcee():
    print('@@@@@@@ Mcor_emcee @@@@@@@')

    for test_key in firstderivs_functions:

        if test_key:
            print('---',test_key)
            if test_key == 'waterfall':
                names = ['Mcorr',
                         'matter', 'radiation', 'a_ombar', 'b_ombar', 'c_ombar',
                         'v_in', 'w_in', 'x_in', 'y_in', 'z_in']
                values = np.array([-19.3,
                                   0.3, 0.025, 0.1, 0.1, 0.1,
                                   0.0, 0.0, 0.0, 0.0, 0.0])
            elif test_key == 'exotic':
                names = ['Mcorr', 'matter', 'radiation', 'gamma', 'zeta']
                values = np.array([-19.3, 0.3, 0.025, 0.0, 0.0])
            elif test_key == 'LCDM':
                names = ['Mcorr', 'matter']
                values = np.array([-19.3, 0.3])
            else:
                names = ['Mcorr', 'matter','gamma']
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

#emcee()

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

