#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 21:45:40 2018

@author: BallBlueMeercat
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
import os.path
import results
import tools
import datasim
import stats


# Number of datapoints to be simulated and number of emcee steps.
npoints, nsteps = 1000, 10000
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

#    Plots for a model.
#datasim.magn({'m':0.3, 'gamma':0.2}, zpicks, 'heaviside_late_int', plot_key=True)
#    Plots for one model with multiple gammas.
#datasim.model_comparison({'m':0.3}, zpicks, 'heaviside_late_int', [-1.45, 0, 0.2])
#   Plots of various models on the same plot. 
#   To check if models overlap, put the two of interest last.
#firstderivs_list = ['heaviside_late_int', 'late_int', 'late_intxde']
#datasim.model_comparison({'m':0.3, 'gamma':-2}, zpicks, firstderivs_list)

firstderivs_functions = [
            'exotic'
#            'late_intxde'
#            ,'heaviside_late_int'
#            ,'late_int'     # limited prior
#            ,'expgamma'     # limited prior
#            ,'txgamma'
#            ,'zxgamma'      
#            ,'gamma_over_z'
#            ,'zxxgamma'     # nan emcee, limited prior
#            ,'gammaxxz'    # nan emcee
#            ,'rdecay_m'     # nan field
#            ,'rdecay_de'
#            ,'rdecay_mxde'  # nan field
#            ,'rdecay'       # limited prior               
#            ,'interacting'  # nan field, limited prior
            ,'LCDM'         # limited prior
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
        datasim.magn(test_params, data_dict, test_key, plot_key=True)
        
    return

#all_modelcheck()


#def quickemcee():
#    print('@@@@@@@ quickemcee @@@@@@@')
##    mag = datasim.noisy_mag(zpicks, mu, sigma, data_params, data_key)
##    dataname = 'Amanullah_sorted1'
#    mag, zpicks = results.load('./data', dataname)
#    
##    pantheon = pd.read_csv('./data/lcparam_full_long.txt', sep=" ")
##    
##    # Reading each txt file column of interest as numpy.ndarray
##    mag = pantheon.mb.values
##    zpicks = pantheon.zhel.values
##    # Stacking them together and sorting by accending redshift.
##    data = np.stack((mag,zpicks), axis=0)
##    data.sort(axis=-1)
##    mag = data[0]
##    zpicks = data[1]
##    zpicks = zpicks.tolist()
#    
#    data_dict = {'mag':mag, 'zpicks':zpicks}
#    test_params = {'m':0.3, 'gamma':0}
#
#    for test_key in firstderivs_functions:
#        
#        # Creating a folder for saving output.
#        save_path = './quick_emcee/'+str(int(time.time()))+'_'+test_key
#        if not os.path.exists(save_path):
#            os.makedirs(save_path)    
#      
#        # Script timer.
#        timet0 = time.time()            
#        
#        # emcee parameter search.
#        propert, sampler = stats.stats(test_params, data_dict, sigma, 
#                                       nsteps, save_path, test_key)        
#        # Time taken by script. 
#        timet1=time.time()
#        tools.timer('script', timet0, timet1)
#        
#        # Saving sampler to directory.
#        results.save(save_path, 'sampler', sampler)
#
#        print('Model being tested:', test_key)
#        print('Data:',dataname)
#
#    return
#
#quickemcee()
    
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
    test_params = {'m':0.3,'M':-19.3,'alpha':0,'beta':0,'gamma':0,'zeta':0}

    
#    mag, zpicks = results.load('./data', dataname)
#    data_dict = {'mag':mag,'zpicks':zpicks}
#    test_params = {'m':0.3, 'M':-19.3, 'gamma':0}

    for test_key in firstderivs_functions:
        
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

        print('Model being tested:', test_key)
        print('Data:',dataname)

    return

Mcor_emcee()


def errorvsdatasize():
    
    data_key = 'LCDM'
    test_key = 'late_int'
    
    data_params = {'m':0.3, 'gamma':0}
    test_params = {'m':0.3, 'gamma':0}
    
    # Script timer.
    timet0 = time.time()
    
    sigma = 0.02
    sigma_max = 0.03
    sigma_step = 0.05
    npoints_min = 1000
    npoints_max = 1100
    npoints_step = 3000
    
    # How many iterations have I signed up for?
    tools.runcount(sigma, sigma_max, sigma_step,
              npoints_min, npoints_max, npoints_step)
    
    decision = input('Happy with the number of iterations? (enter=yes) ')
    if  decision:
        return
    
    # Folder for saving output.
    directory = str(int(time.time()))
    # Relative path of output folder.
    save_path = './emcee_results/'+directory
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    
    run = 0
    
    m_sd_l = []
    m_mean_l = []
    m_vc_l = []    

    g_sd_l = []
    g_mean_l = []
    g_vc_l = [] 
    
    sigma_l = []
    npoints_l = []
    sampler_l = []
    
    while sigma < sigma_max:

        npoints = npoints_min 
        
        # Data to be used:
        mag = datasim.noisy_mag(mu, sigma, npoints, data_params, data_key)
        
        while npoints < npoints_max:
            print('_____________________ run number',run)
          
            propert, sampler = stats.stats(test_params, zpicks, mag, 
                                           sigma, nsteps, save_path, test_key)
            
            m_sd = propert.get('m_sd',0)
            m_mean = propert.get('m_mean', 0)
            m_vc = m_sd/m_mean * 100
            m_vc_l.append(m_vc)
            m_sd_l.append(m_sd)
            m_mean_l.append(m_mean)
            
            g_sd = propert.get('gamma_sd', 0)
            g_mean = propert.get('gamma_mean', 0)
            g_vc = g_sd/g_mean * 100
            g_vc_l.append(g_vc)
            g_sd_l.append(g_sd)
            g_mean_l.append(g_mean)                        
            
            sigma_l.append(sigma)
            npoints_l.append(npoints)
            sampler_l.append(sampler)
            
            npoints += npoints_step
            run += 1
        
        sigma += sigma_step
        
    # Saving plots to run directory.
    # m
    plt.figure()
    plt.xlabel('size of dataset')
    plt.ylabel('standard deviation of marginalised m distribution')
    plt.title('sd of m vs size of dataset, sd of noise = %s'%(sigma))
    plt.scatter(npoints_l, m_sd_l, c='m')        
    stamp = str(int(time.time()))
    filename = str(stamp)+'_sd_of_m_.pdf'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)
    plt.show()
    
    plt.figure()
    plt.xlabel('size of dataset')
    plt.ylabel('mean of marginalised m distribution')
    plt.title('mean of m vs size of dataset, sd of noise = %s'%(sigma))
    plt.scatter(npoints_l, m_mean_l, c='c')        
    stamp = str(int(time.time()))
    filename = str(stamp)+'_mean_of_m_.pdf'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)
    plt.show()
    
    plt.figure()
    plt.xlabel('size of dataset')
    plt.ylabel('variance coefficient in %')
    plt.title('sd/mean of m vs size of dataset, sd of noise = %s'%(sigma))
    plt.scatter(npoints_l, m_vc_l, c='coral')        
    stamp = str(int(time.time()))
    filename = str(stamp)+'_cv_of_m_.pdf'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)
    plt.show()

    # gamma
    plt.figure()
    plt.xlabel('size of dataset')
    plt.ylabel('standard deviation of marginalised gamma distribution')
    plt.title('sd of gamma vs size of dataset, sd of noise = %s'%(sigma))
    plt.scatter(npoints_l, g_sd_l, c='m')        
    stamp = str(int(time.time()))
    filename = str(stamp)+'_sd_of_g_.pdf'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)
    plt.show()
    
    plt.figure()
    plt.xlabel('size of dataset')
    plt.ylabel('mean of marginalised gamma distribution')
    plt.title('mean of gamma vs size of dataset, sd of noise = %s'%(sigma))
    plt.scatter(npoints_l, g_mean_l, c='c')        
    stamp = str(int(time.time()))
    filename = str(stamp)+'_mean_of_g_.pdf'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)
    plt.show()
    
#    plt.figure()
#    plt.xlabel('size of dataset')
#    plt.ylabel('variance coefficient in %')
#    plt.title('sd/mean of gamma vs size of dataset, sd of noise = %s'%(sigma))
#    plt.scatter(npoints_l, g_vc_l, c='coral')        
#    plt.stamp = str(int(time.time()))
#    plt.filename = str(stamp)+'_cv_of_g_.pdf'
#    plt.filename = os.path.join(save_path, filename)
#    plt.savefig(filename)
#    plt.show()
        
    # Saving results to directory.
    results.save(save_path, 'm_vc', m_vc_l)
    results.save(save_path, 'm_sd', m_sd_l)
    results.save(save_path, 'm_mean', m_mean_l)

    results.save(save_path, 'g_vc', g_vc_l)
    results.save(save_path, 'g_sd', g_sd_l)
    results.save(save_path, 'g_mean', g_mean_l)
    
    results.save(save_path, 'sigma', sigma_l)
    results.save(save_path, 'npoints', npoints_l)
    results.save(save_path, 'sampler', sampler_l)
    
    print('directory:',directory)
    
    # Time taken by evaluator. 
    timet1=time.time()
    tools.timer('evaluator', timet0, timet1)
    
    return

#errorvsdatasize()