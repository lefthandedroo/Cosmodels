#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 19:05:26 2018

@author: BallBlueMeercat
"""
import numpy as np 
import matplotlib.pyplot as plt
import time
import os
import zodesolve
import tools
import datasim
import results
#datasim

def magn(params, data, firstderivs_key, plot_key=False):
    """
    Finding matter density m, interaction gamma.
    
    Takes in:
            params = dictionary with true parameters;
            zpicks = list of redshifts to integrate over, in accending order;
            firstderivs_key = string, indicates which firstderivs to integrate;
            plot_key = Boolean, to plot or not to plot model figures;
    Returns:
        mag = np.ndarray of apparent mag corresponding to input redshits.
    """
#    print('@@@ magn has been called')
    if firstderivs_key == 'LCDM':
        params['gamma'] = 0
        del params['gamma']
    
    zpicks = data['zpicks']
    
    # Absolute brightness of supernovae.
    M = -19
    
    dlpc, plot_var = zodesolve.zodesolve(params, zpicks, firstderivs_key)
    
    # Calculating apparent magnitudes of supernovae at the simulated
    # luminosity distances using the distance modulus formula.
    mag = 5 * np.log10(dlpc/10) + M

    if plot_key:
        # Checking evolution of the model.
        import plots
        plots.modelcheck(mag, zpicks, plot_var, firstderivs_key)
        
    return mag

#def magn(params, data, firstderivs_key, plot_key=False):
#    """
#    Finding matter density m, alpha, beta, interaction gamma.
#    Takes in:
#            params = dictionary with true parameters;
#            zpicks = list of redshifts to integrate over, in accending order;
#            firstderivs_key = string, indicates which firstderivs to integrate;
#            plot_key = Boolean, to plot or not to plot model figures;
#    Returns:
#        mag = np.ndarray of apparent mag corresponding to input redshits.
#    """
##    print('@@@ magn has been called')
#    if firstderivs_key == 'LCDM':
#        params['gamma'] = 0
#        del params['gamma']
#    
#    zpicks = data['zpicks']
#    x1 = data['x1']
#    colour = data['colour']
#    
#    # Absolute brightness of supernovae.
#    M = params['M']
#    alpha = params['alpha']
#    beta = params['beta']
#    
#    dlpc, plot_var = zodesolve.zodesolve(params, zpicks, firstderivs_key)
#    
#    # Calculating apparent magnitudes of supernovae at the simulated
#    # luminosity distances using the distance modulus formula.
#    mag = 5 * np.log10(dlpc/10) + M - alpha*x1 +beta*colour
#
#    if plot_key:
#        # Checking evolution of the model.
#        import plots
#        plots.modelcheck(mag, zpicks, plot_var, firstderivs_key)
#        
#    return mag
    
# Slow mag calculation
#    # Calculating apparent magnitudes of supernovae at the simulated
#    # luminosity distances using the distance modulus formula.
#    mag = []   
#    for i in range(len(dlpc)):
#        if dlpc[i] == 0:
#            magnitude = M
#        else:
#            # magnitude from the distance modulus formula
#            magnitude = 5 * math.log10(dlpc[i]/10) + M
#        mag.append(magnitude)
    
#evaluator

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