3#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:47:44 2018

@author: BallBlueMeercat
"""
import random
import math
import numpy as np
import zodesolve
import tools

def redshift_picks(zmin, zmax, n):
    """
    Takes in:
        zmin = integer lowest redshift;
        zmax = integer highest redshift;
        n = integer number of redshifts to be generated.
    Returns:
        zpicks = list of randomly selected redshifts between zmin and zmax.
    """
#    print('-zpicks has been called')
    zinterval = (zmax - zmin) / (n*2)
    z_opts = tools.flist(zmin, zmax, zinterval)
    zpicks = random.sample(z_opts, n)
    zpicks = sorted(zpicks)
    return zpicks

def magn(params, zpicks, firstderivs_key, plot_key=False):
    """
    Takes in:
            params = dictionary with true parameters;
            zpicks = list of redshifts to integrate over;
            firstderivs_key = string, indicates which firstderivs to integrate;
            plot_key = Boolean, to plot or not to plot model figures;
    Returns:
        mag = np.ndarray of apparent mag corresponding to input redshits.
    """
#    print('@@@ zmsim has been called')
    
    # Absolute brightness of supernovae.
    M = -19
    
    if not sorted(zpicks) == zpicks:
        zpicks.sort()
        print('sorted to accending in zmsim')    
    
    if not isinstance(zpicks, (list,)):
        zpicks = zpicks.tolist()
        print('converted to list in zmsim')
    
    dlpc, plot_var = zodesolve.zodesolve(params, zpicks, firstderivs_key)
    
    # Calculating apparent magnitudes of supernovae at the simulated
    # luminosity distances using the distance modulus formula.
    mag = 5 * np.log10(dlpc/10) + M

    if plot_key:
        # Checking evolution of the model.
        import plots
        plots.modelcheck(mag, zpicks, plot_var, firstderivs_key)
        
    return mag

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
    
#def magn(params, zpicks, firstderivs_key, plot_key=False, gamma_list=False):
#    """
#    Takes in:
#            gamma = interaction constant;
#            m = e_m(t)/ec(t0) at t=t0;
#            de = e_de(t)/ec(t0) at t=t0;
#            zpicks = list of z to match the interpolated dlmpc to;
#            gamma_list = list of gamma to subplot on model figures
#    Returns:
#        mag = list of n apparent magnitudes mag corresponding to given redshits.
#    """
##    print('@@@ zmsim has been called')
#    
#    # Absolute brightness of supernovae.
#    M = -19
#    
#    if not sorted(zpicks) == zpicks:
#        zpicks.sort()
#        print('sorted to accending in zmsim')    
#    
#    if not isinstance(zpicks, (list,)):
#        zpicks = zpicks.tolist()
#        print('converted to list in zmsim')
#        
#    if gamma_list:
#        plot_var_dict = {}
#
#        j = 1        
#        for gamma in gamma_list:
#            params['g_true'] = gamma
#            dlpc, plot_var = zodesolve.zodesolve(params, zpicks, firstderivs_key)
#        
#            # Calculating apparent magnitudes of supernovae at the simulated
#            # luminosity distances using the distance modulus formula.
#            mag = 5 * np.log10(dlpc/10) + M
#            
#            plot_var_dict['plot_var_'+str(j)] = plot_var
#            print('plot_var_'+str(j))
#            print(j)
#            plot_var_dict['mag_'+str(j)] = mag
#            
#            j+=1
#        
#        if plot_key:
#            # Checking evolution of the model.
#            import plots
#            plots.modelcheck(mag, zpicks, plot_var, firstderivs_key, plot_var_dict)
#            
#    else:
#        dlpc, plot_var_1 = zodesolve.zodesolve(params, zpicks, firstderivs_key)
#        
#        # Calculating apparent magnitudes of supernovae at the simulated
#        # luminosity distances using the distance modulus formula.
#        mag = 5 * np.log10(dlpc/10) + M
#        
#        if plot_key:
#            # Checking evolution of the model.
#            import plots
#            plots.modelcheck(mag, zpicks, plot_var_1, firstderivs_key)
#        
#    return mag


def gnoise(mag, mu, sigma):
    """
   Returns:
       mag = mag, each point offset by unique Gaussian noise;
       noise = Gaussian noise.
    """
#    print('               -gnoise has been called')
    n = len(mag)
    noise = np.random.normal(mu,sigma,n)
    mag = mag + noise
    
#    import matplotlib.pyplot as pl
#    from pylab import figure
#    figure()
#    pl.title('Noise distribution')
#    pl.hist(noise, 100)
#    pl.show()
    return mag


def noisy_mag(zpicks, mu, sigma, params, firstderivs_key, plot_key=False):
    
    model = magn(params, zpicks, firstderivs_key, plot_key)
    model = np.asarray(model)
    mag = gnoise(model, mu, sigma)
    
    return mag


def makensavemagnz(m_true, g_true, mu, sigma, zpicks, data_key, filename):
    '''
    Takes in:
    
        Parameters used to simulate magnitude:  
    m_true                  = e_m(t)/e_crit(t0) at t=t0;
    de_true = 1 - m_true    = e_de(t)/e_crit(t0) at t=t0;
    g_true = interaction term, rate at which DE decays into matter.
    
        Statistical parameteres of gaussian noise added to data:
    mu =  mean;
    sigma = standard deviation;
    
    npoints = how many mag and z to generate.
    
        Model type:
    data_key = string, key for dictionary of interaction modes in firstderivs
    Options: 'Hdecay', 'rdecay', 'rdecay_de', 'rdecay_m', 'interacting', 'LCDM'
    Length of parameters has to correspond to the model being tested.
    
    filename = string, name of file data is saved to.
    
    Returns:
        Nothing. Generates redshifts and corresponding magnitudes (according 
        to the model specified by data_key) offset by Gaussian noise, 
        saves them into a binary file called filename in the working directory.
    '''
    
    if data_key == 'LCDM':
        data_params = {'m':m_true}
    else:
        data_params = {'m':m_true, 'gamma':g_true}
    
    mag = noisy_mag(zpicks, mu, sigma, data_params, data_key)
    
    output = mag, zpicks
    
    # Relative path of output folder.
    save_path = './data/'+filename
        
    import pickle
    pickle.dump(output, open(save_path, 'wb'))
    
    return





