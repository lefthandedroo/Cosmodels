#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 21:45:40 2018

@author: BallBlueMeercat
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('default') # has to be switched on to set figure size
mpl.style.use('fivethirtyeight')
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['grid.color'] = 'white'

import numpy as np
import time
import pickle
import os.path
import results
import tools
import datasim
import stats
print('- - - - - - - - evaluator')
timed = False
plot = True # emcee plots

# Number of emcee steps.
nsteps = 10000

# Statistical parameteres of noise: mean, standard deviation.
mu, sigma = 0.0, 0.001 # sigma != 0

data_dic = {}

dataname = 'pantheon'
#dataname = 'LCDM_to_1089'
#dataname = 'LCDM_to_2.26'
#dataname = 'specific_z'
if dataname == 'pantheon':
    import pandas as pd
    print('-----Using pantheon')
    # Pantheon data:
    pantheon = pd.read_csv('./data/lcparam_full_long.txt', sep=" ")
    pantheon.set_index('name', inplace=True)
    pantheon.sort_values('zhel', inplace=True)
    mag = pantheon.mb.values
    sigma = pantheon.dmb.values
    zpicks = pantheon.zhel.values
#    plt.figure()
#    plt.hist(sigma, 70, facecolor="C{}".format(5), label='pantheon errors')
#    plt.xlabel(r'$\delta$Mag')
#    plt.legend()
#    plt.show()
    data_dic['mag'] = mag
elif dataname == 'LCDM_to_1089':
    print('-----Generating z up to z=1089')
    # Generating artificial redshifts.
    zpicks = np.sort(np.random.uniform(low=0.01012, high=1088, size=(10000,)))
    zpicks[-1] = 1089
    # Generating LCDM mag and da.
    names, values = ['Mcorr', 'matter'], np.array([-19.3, 0.3])
    mag, da = datasim.magn(names, values, {'zpicks':zpicks}, 'LCDM', plot_key=False)
    # Adding noise to LCDM mag.
    nmag = datasim.gnoise(mag, mu, sigma)
    data_dic['mag'] = nmag
elif dataname == 'LCDM_to_2.26':
    print('-----Generating z up to z=2.26')
    # Generating artificial redshifts.
    zpicks = np.sort(np.random.uniform(low=0.01012, high=2.26, size=(1048,)))
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

data_dic['data_zpicks'] = zpicks

#zpicks = np.sort(np.random.uniform(low=0.01012, high=zpicks[-1], size=(1000,)))
#zpicks[-1] = zpicks[-1]
data_dic['zpicks'] = zpicks

model_name = ['rdecay']
p0 = ['Mcorr', 'm_ombar', 'gamma'], np.array([-19.3, 0.3, -1])
p1 = ['Mcorr', 'm_ombar', 'gamma'], np.array([-19.3, 0.3, 1])
#datasim.model_comparison([p0, p1], data_dic, model_name*2, plot_key=True)








# Plot param evolutions for multiple models on the same axis.
model_name = ['interacting', 'interacting', 'rdecay_m', 'zxgamma', 'expgamma', 'heaviside_late_int', 'late_intxde']#, 'gamma_over_z']
names = ['Mcorr','matter','gamma']

p2 = names, np.array([-19.3, 0.28, 0])      # LCDM
p3 = names, np.array([-19.31, 0.2, -0.52])  # interacting
p4 = names, np.array([-19.3, 0.22, -0.65])  # rdecay_m
p5 = names, np.array([-19.3, 0.28, -0.34])  # zxgamma
p6 = names, np.array([-19.3, 0.24, 0.56])   # expgamma
p7 = names, np.array([-19.31, 0.2, -0.52])  # heaviside_late_int
p8 = names, np.array([-19.3, 0.26, -0.57])  # late_intxde
#p9 = names, np.array([-19.31, 0.15, -0.23]) # gamma_over_z
datasim.model_comparison([p2, p3, p4, p5, p6, p7, p8], data_dic, model_name, plot_key=True)











#model_name = ['exotic']
#names = ['Mcorr','matter','radiation', 'gamma','zeta']
#p2 = names, np.array([-19.3, 0.3, 0.025, 1, 1])
#p3 = names, np.array([-19.3, 0.3, 0.025, 0.7, 0.7])
#p4 = names, np.array([-19.3, 0.3, 0.025, 0.4, 0.4])
#p5 = names, np.array([-19.3, 0.3, 0.025, -0.01, -0.01])
#datasim.model_comparison([p2, p3, p4, p5], data_dic, model_name*4, plot_key=True)

#p1 = ['Mcorr','m_ombar', 'r_ombar', 'a_ombar', 'b_ombar', 'c_ombar',
#     'd_ombar', 'e_ombar', 'f_ombar', 'g_ombar', 'h_ombar',
#     'i_ombar','a_in', 'b_in', 'c_in', 'd_in', 'e_in', 'f_in',
#     'g_in', 'h_in', 'i_in', 'j_in', 'k_in'], np.array([-19.3,0.3, 0.025, 0, 0,
#    0, 0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
#datasim.model_comparison([p1], data_dic, ['rainbow'], plot_key=True)

firstderivs_functions = [None
            ,'rainbow'
            ,'niagara'
            ,'kanangra'
            ,'waterfall'
            ,'stepfall'
            ,'exotic'
            ,'late_intxde'
            ,'heaviside_late_int'
            ,'heaviside_sudden'
            ,'late_int'
            ,'expgamma'
            ,'txgamma'
            ,'zxgamma'
            ,'gamma_over_z'
            ,'zxxgamma'        # gamma forced positive in firstderivs
            ,'gammaxxz'        # gamma forced positive in firstderivs
            ,'rdecay_m'
            ,'rdecay_de'
            ,'rdecay_mxde'
            ,'rdecay'
            ,'interacting'
            ,'LCDM'
            ,'rLCDM'
            ]

if timed:
    import cProfile, pstats, io
    pr = cProfile.Profile()
    pr.enable()

def modelcheck():
    for test_key in firstderivs_functions:
        if test_key:
            print('---',test_key)
            names, values, int_in = tools.names_values(test_key)
            datasim.magn(names, values, data_dic, test_key, plot=True)
    return

#modelcheck()

def emcee():
    print('@@@@@@@ emcee @@@@@@@')

    for test_key in firstderivs_functions:
        if test_key:
            print('---',test_key)
            names, values, int_in = tools.names_values(test_key)

            # Creating a folder for saving output.

            save_path = './results_emcee/'+str(int(time.time()))+'_'+test_key
            if not os.path.exists(save_path):
                    os.makedirs(save_path)

            # Script timer.
            timet0 = time.time()

            # emcee parameter search.
            propert, sampler = stats.stats(names, values, data_dic,
                                           sigma, nsteps, save_path,
                                           test_key, int_in, plot=plot)
            # Time taken by script.
            timet1=time.time()
            tools.timer('script', timet0, timet1)

            # Saving sampler to directory.
            results.save(save_path, 'sampler', sampler)

            print('Data:',dataname)

    return

#emcee()

if timed:
    pr.disable()
    s = io.StringIO()
    sortby = 'tottime'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print (s.getvalue())