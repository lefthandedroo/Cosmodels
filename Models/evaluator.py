#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 21:45:40 2018

@author: BallBlueMeercat
"""
import numpy as np
import time
import pickle
import os.path
import results
import tools
import datasim
import stats


# Number of emcee steps.
nsteps = 10000

# Statistical parameteres of noise: mean, standard deviation.
mu, sigma = 0.0, 0.007 # sigma != 0

dataname = 'pantheon'
if dataname == 'pantheon':
    import pandas as pd
    print('-----Using pantheon')
    # Pantheon data:
    pantheon = pd.read_csv('./data/lcparam_full_long.txt', sep=" ")
    pantheon.set_index('name', inplace=True)
    pantheon.sort_values('zhel', inplace=True)
    mag = pantheon.mb.values
    zpicks = pantheon.zhel.values
elif dataname == 'generated synth':
    ######### Artificial LCDM data: ##########
    # Generating redshifts.
    zpicks = np.sort(np.random.uniform(low=0.0001, high=1088, size=(10000,)))
    zpicks[-1] = 1089
elif dataname == 'pre-made zpicks':
    # Extracting pre-made redshifts z=0 to z=1089.
    try:
        with open('data/zpicks_1000_1089.p','rb') as rfp: zpicks = pickle.load(rfp)
    except:
        print("zpicks_1000_1089.p didnt't open")

data_dic = {'zpicks':zpicks}

# Generating LCDM mag and da.
names, values = ['Mcorr', 'matter'], np.array([-19.3, 0.3])
mag, da = datasim.magn(names, values, data_dic, 'LCDM', plot_key=False)

# Adding noise to LCDM mag.
nmag = datasim.gnoise(mag, mu, sigma)
dataname = 'noisy LCDM'
data_dic['mag'] = nmag

## Plot param evolutions for multiple models on the same axis.
#p1 = ['Mcorr', 'm_ombar'], np.array([-19.3, 0.3])
#p2 = ['Mcorr', 'm_ombar', 'r_ombar'], np.array([-19.3, 0.3, 0.025])
#datasim.model_comparison([p1, p2], zpicks, ['LCDM', 'rLCDM'], plot_key=True)


firstderivs_functions = [None
#            ,'stepfall'
#            ,'waterfall'
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
#            ,'LCDM'
            ,'rLCDM'
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
            elif test_key == 'rLCDM':
                names = ['Mcorr', 'm_ombar', 'r_ombar']
                values = np.array([-19.3, 0.3, 0.025])
            elif test_key == 'LCDM':
                names = ['Mcorr', 'm_ombar']
                values = np.array([-19.3, 0.3])
            else:
                names = ['Mcorr', 'm_ombar','gamma']
                values = np.array([-19.3, 0.3, 0.0])

            # Making sure number of parameters matches number of names given:
            assert len(names) == len(values), "len(names) != len(values)"
            datasim.magn(names, values, data_dic, test_key, plot_key=True)
    return

modelcheck()

def emcee():
    print('@@@@@@@ emcee @@@@@@@')

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
            assert len(names) == len(values), "len(names) != len(values)"

            # Creating a folder for saving output.
            save_path = './results_emcee/'+str(int(time.time()))+'_'+test_key
            if not os.path.exists(save_path):
                os.makedirs(save_path)

            # Script timer.
            timet0 = time.time()

            # emcee parameter search.
            propert, sampler = stats.stats(names, values, data_dic,
                                           sigma, nsteps, save_path,
                                           test_key, plot=1)
            # Time taken by script.
            timet1=time.time()
            tools.timer('script', timet0, timet1)

            # Saving sampler to directory.
            results.save(save_path, 'sampler', sampler)

            print('Data:',dataname)

    return

#emcee()