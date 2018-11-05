#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 16:02:10 2018

@author: BallBlueMeercat
"""
import matplotlib.pyplot as plt
from emcee import EnsembleSampler
import numpy as np
import time
import os.path

import datasim
import tools
import ln
import plots

def stats(test_params, data_dict, sigma, nsteps, 
          save_path, firstderivs_key, plot=False):
    """
    Takes in:
        test_params = list of dictionaries {string:value} of names and 
        starting values of parameters to be emcee fitted:
            {'matter':int/float} = e_m(t)/ec(t0) at t=t0;
            {'Mcorr':int/float} = corrected absolute mag M;
            {'gamma':int/float} = interaction term;
            {'zeta':int/float} = interaction term;
            ... (more)
        data_dict = dictionary of parameters from data:
            {'colour': numpy.ndarray} = SN colour;
            {'x1': numpy.ndarray} = SN stretch correction as;
            {'zpicks':list} = redshifts sorted in accending order;
            {'mag':list} = apparent magnitudes;
        sigma = standard deviation of error on the data;
        nsteps = int, steps to be taken by each emcee walker;
        save_path = string, directory for saving output;
        firstderivs_key = string, name of model to use for model mag.
    Returns:
        propert, sampler
    """
    
    zpicks = data_dict.get('zpicks',0)
    mag = data_dict.get('mag',0)
    
    if firstderivs_key == 'rainbow' or firstderivs_key == 'waterfall':
        test_params = test_params
    elif firstderivs_key == 'exotic':
        test_params = test_params[:4]
    elif firstderivs_key == 'LCDM':
        test_params = test_params[:2]
    else:
        test_params = test_params[:3]
    
    # emcee parameters:
    ndim = len(test_params)
    nwalkers = int(ndim * 2)
    
    # Initializing walkers.
    pos = np.zeros((1, ndim)).flatten()
    for i in range(ndim):
        dic = test_params[i]
        for key in dic:
            pos[i] = dic[key]
    pos = [pos + 0.001*np.random.randn(ndim) for i in range(nwalkers)]
    
    # Are walkers starting outside of prior?
    for i in range(nwalkers):
        theta = pos[i]
        lp = ln.lnprior(theta, firstderivs_key)
        if not np.isfinite(lp):
            print('~~~~~~~pos[%s] (outside of prior) = %s ~~~~~~~'%(i, theta))
        
    # Sampler setup.
    times0 = time.time()    # starting sampler timer
    sampler = EnsembleSampler(nwalkers, ndim, ln.lnprob, 
                              args=(data_dict, sigma, firstderivs_key, test_params))
    
    # Burnin.
    burnin = int(nsteps/4)  # steps to discard
    print('_____ burnin start')
    timeb0 = time.time()    # starting burnin timer
    pos, prob, state = sampler.run_mcmc(pos, burnin)
    timeb1=time.time()      # stopping burnin timer
    print('_____ burnin end')
    sampler.reset()
    
    # Starting sampler after burnin.
    print('_____ sampler start')
    sampler.run_mcmc(pos, nsteps)
    print('_____ sampler end')
    times1=time.time()      # stopping sampler timer
    
    # Walker steps is lnprob = sampler.flatlnprobability
    # Index of parameter values with highest post prob found by emcee.
    bi = np.argmax(sampler.flatlnprobability)
    
    # Extracting results:
    thetabest = np.zeros(ndim)
    true = []
    truth_names = []
    propert = {}
    propert['trace'] = sampler.chain[:, burnin:, :].reshape(-1, ndim) #trace
    
    colours = ['light red', 'berry', 'coral', 'amber', 
               'apple', 'aquamarine', 'raspberry', 'green blue', 'deep blue',
               'emerald', 'blue violet', 'dark violet', 'yellow orange']    
    
    best_params = test_params
    
    for i in range(ndim):
        
        for key in test_params[i]:
            truth_names.append(key)
            param_initial = key[0]
            best = sampler.flatchain[bi,i]
            # Input parameter.
            param_true = test_params[i].get(key, 0)
            true.append(param_true)
            # Output parameter.
            output = sampler.flatchain[:,i]
            # Standard deviation and mean of the emcee found distribution.
            propert[param_initial+'_sd'] = np.std(output)
            propert[param_initial+'_mean'] = np.mean(output)
            propert[param_initial] = sampler.flatchain[bi,i]
            
            if plot:
                plots.stat_emcee(colours[i], output, param_true, key, 
                           sampler.flatlnprobability, zpicks, mag, sigma,
                           nsteps, nwalkers, save_path, firstderivs_key)
            
            thetabest[i] = best
            best_params[i][key] = best
            
    # Checking if best found parameters are within prior.
    lp = ln.lnprior(thetabest, firstderivs_key)
    if not np.isfinite(lp):
        print('~~~~~~~thetabest outside of prior at magbest~~~~~~~')

    if plot:
        # Plot of data mag and redshifts, overlayed with
        # mag simulated using emcee best parameters and data redshifts.
        magbest = datasim.magn(best_params, data_dict, firstderivs_key)
        plt.figure()
        plt.title('model: '+firstderivs_key
                  +'\n Evolution of magnitude with redshift \n nsteps: '
              +str(nsteps)+', noise: '+str(sigma)+', npoints: '+str(len(zpicks)))
        data = plt.errorbar(zpicks, mag, yerr=sigma, fmt='.', alpha=0.3)
        best_fit = plt.scatter(zpicks, magbest, lw='1', c='xkcd:tomato')
        plt.ylabel('magnitude')
        plt.xlabel('z')
        plt.legend([data, best_fit], ['LCDM', firstderivs_key])
        stamp = str(int(time.time()))
        filename = str(stamp)+'____magz__nsteps_'+str(nsteps)+'_nwalkers_' \
        +str(nwalkers)+'_noise_'+str(sigma)+'_numpoints_'+str(len(zpicks))+'.png'
        filename = os.path.join(save_path, filename)
        plt.savefig(filename)
    
        # Corner plot (walkers' walk + histogram).
        import corner
        samples = sampler.chain[:, :, :].reshape((-1, ndim))
        corner.corner(samples, labels=truth_names, truths=true)
        plt.show(block=False)
    
    # Results getting printed:
    if bi == 0: 
        print('@@@@@@@@@@@@@@@@@')
        print('best index =',str(bi))
        print('@@@@@@@@@@@@@@@@@')
    print('max likelihood params =',str(best_params))
    print('m.a.f.:', np.mean(sampler.acceptance_fraction))
    print('nsteps:', str(nsteps))
    print('sigma:', str(sigma))
    print('npoints:', str(len(zpicks)))
    print('model:', firstderivs_key)
    
    tools.timer('burnin', timeb0, timeb1)
    tools.timer('sampler', times0, times1)
    
    return propert, sampler