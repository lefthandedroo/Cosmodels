#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 21:45:40 2018

@author: BallBlueMeercat
"""
import pandas as pd
import numpy as np
import time
import os.path
import results
import tools
import datasim
import stats


# Number of datapoints to be simulated and number of emcee steps.
npoints, nsteps = 1000, 1000
zmax = 2

# Statistical parameteres of noise:
mu = 0.0            # mean
sigma = 0.07      # standard deviation

dataname = 'mag_z_LCDM_1000_sigma_'+str(sigma)

#    Making data (mag and z)
#dataname = 'mag_z_'+data_key+'_'+str(npoints)+'_sigma_'+str(sigma)
#datasim.makensavemagnz(m_true, g_true, mu, sigma, zpicks, data_key, dataname)

#    Making redshifts to use in this script.
zpicks = datasim.redshift_picks(0.005, zmax, npoints)

#    Plots for one model with multiple gamma or zeta.
#datasim.model_comparison({'m':0.3,'gamma':0.7}, zpicks, 'expgamma', gamma_list = [-0.3, 0, 0.6])

#   Plots of various models on the same plot. 
#   To check if models overlap, put the two of interest last.
firstderivs_list = ['expgamma', 'exotic', 'rdecay']
datasim.model_comparison({'m':0.3, 'gamma':-0.1}, zpicks, firstderivs_list)

firstderivs_functions = [None
            ,'exotic'
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
             ]

def all_modelcheck():
    print('@@@@@@@ all_modelcheck @@@@@@@')
    
    dataname = './data/lcparam_full_long.txt'
    pantheon = pd.read_csv(dataname, sep=" ")
    
    # Reading each txt file column of interest as numpy.ndarray
    mag = pantheon.mb.values
    x1 = pantheon.x1.values
    colour = pantheon.color.values
    zpicks = pantheon.zhel.values
    
    # Stacking them together and sorting by accending redshift.
    data = np.stack((mag,x1,colour,zpicks), axis=0)
    data.sort(axis=-1)
    
    mag = data[0]
    x1 = data[1]
    colour = data[2]
    zpicks = data[3]
    zpicks = zpicks.tolist()
    data_dict = {'mag':mag, 'x1':x1, 'colour':colour, 'zpicks':zpicks}
    test_params = {'m':0.3,'M':-19.3,'alpha':0,'beta':0,'gamma':0,'zeta':0}
    
    for test_key in firstderivs_functions:
        if test_key:
            datasim.magn(test_params, data_dict, test_key, plot_key=True)    
    return

#all_modelcheck()

    
def Mcor_emcee():
    print('@@@@@@@ Mcor_emcee @@@@@@@')
    dataname = './data/lcparam_full_long.txt'
    pantheon = pd.read_csv(dataname, sep=" ")
    
    # Reading each txt file column of interest as numpy.ndarray
    mag = pantheon.mb.values
    x1 = pantheon.x1.values
    colour = pantheon.color.values
    zpicks = pantheon.zhel.values
    
    # Stacking them together and sorting by accending redshift.
    data = np.stack((mag,x1,colour,zpicks), axis=0)
    data.sort(axis=-1)
    
    mag = data[0]
    x1 = data[1]
    colour = data[2]
    zpicks = data[3]
    zpicks = zpicks.tolist()
    data_dict = {'mag':mag, 'x1':x1, 'colour':colour, 'zpicks':zpicks}
    test_params = [{'matter':0.3},{'Mcorr':-19.3},{'alpha':0},{'beta':0},{'gamma':0},{'zeta':0}]

    for test_key in firstderivs_functions:
        if test_key:
            print('---',test_key)
            # Creating a folder for saving output.
            save_path = './quick_emcee/'+str(int(time.time()))+'_'+test_key
            if not os.path.exists(save_path):
                os.makedirs(save_path)    
          
            # Script timer.
            timet0 = time.time()            
            
            # emcee parameter search.
            propert, sampler = stats.stats(test_params, data_dict, sigma, 
                                           nsteps, save_path, test_key)        
            # Time taken by script. 
            timet1=time.time()
            tools.timer('script', timet0, timet1)
            
            # Saving sampler to directory.
            results.save(save_path, 'sampler', sampler)
    
            print('Data:',dataname)

    return

#Mcor_emcee()