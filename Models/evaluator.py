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
nsteps = 100000

# Statistical parameteres of noise: mean, standard deviation.
mu, sigma = 0.0, 0.03 # sigma != 0

data_dic = {}

dataname = 'pantheon'
#dataname = 'LCDM_to_1089'
#dataname = 'specific_z'
if dataname == 'pantheon':
    import pandas as pd
    print('-----Using pantheon')
    # Pantheon data:
    pantheon = pd.read_csv('./data/lcparam_full_long.txt', sep=" ")
    pantheon.set_index('name', inplace=True)
    pantheon.sort_values('zhel', inplace=True)
    mag = pantheon.mb.values
    zpicks = pantheon.zhel.values
    data_dic['mag'] = mag
elif dataname == 'LCDM_to_1089':
    print('-----Generating z up to z=1089')
    # Generating artificial redshifts.
    zpicks = np.sort(np.random.uniform(low=0.0001, high=1088, size=(10000,)))
    zpicks[-1] = 1089
    # Generating LCDM mag and da.
    names, values = ['Mcorr', 'matter'], np.array([-19.3, 0.3])
    mag, da = datasim.magn(names, values, {'zpicks':zpicks}, 'LCDM', plot_key=False)
    # Adding noise to LCDM mag.
    nmag = datasim.gnoise(mag, mu, sigma)
    data_dic['mag'] = nmag
elif dataname == 'pre-made_zpicks':
    print('-----Using pre-made z')
    # Extracting pre-made redshifts z=0 to z=1089.
    try:
        with open('data/zpicks_1000_1089.p','rb') as rfp: zpicks = pickle.load(rfp)
    except:
        print("zpicks_1000_1089.p didnt't open")
elif dataname == 'specific_z':
    zpicks = [1, 1.2, 1.2600444250620284, 1.3, 1.4, 1.5, 1.6]
    print(f'-----Using z={zpicks}')

data_dic['zpicks'] = zpicks


#p0 = ['Mcorr', 'm_ombar', 'gamma'], np.array([-19.3, 0.3, -1])
#p1 = ['Mcorr', 'm_ombar', 'gamma'], np.array([-19.3, 0.3, 1])
#datasim.model_comparison([p0, p1], data_dic, ['rdecay', 'rdecay'], plot_key=True)

## Plot param evolutions for multiple models on the same axis.
#p2 = ['Mcorr', 'm_ombar', 'gamma'], np.array([-19.3, 0.3, 1])
#p3 = ['Mcorr', 'm_ombar', 'gamma'], np.array([-19.3, 0.3, -1])
#p4 = ['Mcorr', 'm_ombar', 'gamma'], np.array([-19.3, 0.3, -3])
#p5 = ['Mcorr', 'm_ombar', 'gamma'], np.array([-19.3, 0.3, -6])
#datasim.model_comparison([p2, p3, p4, p5], data_dic, ['interacting', 'interacting', 'interacting', 'interacting'], plot_key=True)

#p1 = ['Mcorr','m_ombar', 'r_ombar', 'a_ombar', 'b_ombar', 'c_ombar',
#     'd_ombar', 'e_ombar', 'f_ombar', 'g_ombar', 'h_ombar',
#     'i_ombar','a_in', 'b_in', 'c_in', 'd_in', 'e_in', 'f_in',
#     'g_in', 'h_in', 'i_in', 'j_in', 'k_in'], np.array([-19.3,0.3, 0.025, 0, 0,
#    0, 0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
#datasim.model_comparison([p1], data_dic, ['rainbow'], plot_key=True)

firstderivs_functions = [None
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

def modelcheck():
    for test_key in firstderivs_functions:
        if test_key:
            print('---',test_key)
            names, values = tools.names_values(test_key)
            datasim.magn(names, values, data_dic, test_key, plot_key=True)
    return

#modelcheck()

def emcee():
    print('@@@@@@@ emcee @@@@@@@')

    for test_key in firstderivs_functions:

        if test_key:
            print('---',test_key)
            names, values = tools.names_values(test_key)

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

emcee()