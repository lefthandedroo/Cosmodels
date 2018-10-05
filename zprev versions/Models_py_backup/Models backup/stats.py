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
          save_path, firstderivs_key):
    """
    Takes in:
            test_params = dictionary of parameters to be emcee fitted
                'm':int/float = e_m(t)/ec(t0) at t=t0;
                'gamma':int/float = interaction term;
                'zeta':int/float = interaction term;
                'alpha':int/float = SN peak mag correlation parameter;
                'beta' :int/float = SN peak mag correlation parameter;
            data_dict = dictionary of parameters from data
                'colour': numpy.ndarray = SN colour;
                'x1': numpy.ndarray = SN stretch correction as;
                'zpicks':list of redshifts sorted in accending order;
                'mag':list of apparent magnitudes;
            sigma = standard deviation of error on the data;
            nsteps = int, steps to be taken by each emcee walker;
            save_path = string, directory for saving output;
            firstderivs_key = string, name of IVCDM model to use for model mag.
    Returns:
    """
#    print('-stats has been called')
    
    zpicks = data_dict.get('zpicks',0)
    mag = data_dict.get('mag',0)
    
    if firstderivs_key == 'exotic':
        pass
    elif firstderivs_key == 'LCDM':
        test_params['gamma'] = 0
        del test_params['gamma']
        test_params['zeta'] = 0
        del test_params['zeta']
    else:
        test_params['zeta'] = 0
        del test_params['zeta']
    
    # emcee parameters:
    ndim = len(test_params)
    nwalkers = int(ndim * 2)
    
    # Initializing walkers.
    poslist = list(test_params.values())
    pos = []
    for i in poslist:
        pos.append(i)
    startpos = np.array(pos)
    pos = [startpos + 0.001*np.random.randn(ndim) for i in range(nwalkers)]
    
    # Are walkers starting outside of prior?
    i = 0
    while i < nwalkers:
        theta = pos[i]
        lp = ln.lnprior(theta, firstderivs_key)
        if not np.isfinite(lp):
            print('~~~~~~~pos[%s] (outside of prior) = %s ~~~~~~~'%(i, theta))
        i += 1
        
    # Sampler setup.
    times0 = time.time()    # starting sampler timer
    sampler = EnsembleSampler(nwalkers, ndim, ln.lnprob, 
                                    args=(data_dict, sigma, firstderivs_key, ndim))
    
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
    
    # Walker steps.
    lnprob = sampler.flatlnprobability
    # Index of best parameters found by emcee.
    bi = np.argmax(sampler.flatlnprobability) # index with highest post prob 
    
    trace = sampler.chain[:, burnin:, :].reshape(-1, ndim)
    
    # Extracting results:
    thetabest = np.zeros(ndim)
    parambest = {}
    true = []
    propert = {}
    propert['trace'] = trace
    
    colours = ['coral', 'orchid', 'apple', 'orange', 'aquamarine', 'black']
    
    def stat(i, sampler, string, test_params, propert):
        
        best_output = sampler.flatchain[bi,i]
        # Input m = e_m(z)/ec(z=0).
        param_true = test_params.get(string, 0)
        true.append(param_true)
        # Output m.
        output = sampler.flatchain[:,i]
        # Standard deviation and mean of the m distribution.
        propert[string+'_sd'] = np.std(output)
        propert[string+'_mean'] = np.mean(output)
        propert[string] = sampler.flatchain[bi,i]
        
        return best_output, output, param_true, propert
    
    
    for i in range(ndim):
        
        if i == 0:            
            best, output, param_true, propert = stat(i, sampler, 'm', test_params, propert)
            
            plots.stat(colours[i], output, param_true, 'matter', lnprob, zpicks, 
                       mag, sigma, nsteps, nwalkers, save_path, firstderivs_key)
            
            thetabest[i] = best
            parambest['m'] = best
                        
        elif i == 1:            
            best, output, param_true, propert = stat(i, sampler, 'M', test_params, propert)
            
            plots.stat(colours[i], output, param_true, 'Mcorr', lnprob, zpicks, 
                       mag, sigma, nsteps, nwalkers, save_path, firstderivs_key)
            
            thetabest[i] = best
            parambest['M'] = best
                        
        elif i == 2:
            best, output, param_true, propert = stat(i, sampler, 'a', test_params, propert)
            
            plots.stat(colours[i], output, param_true, 'alpha', lnprob, zpicks, 
                       mag, sigma, nsteps, nwalkers, save_path, firstderivs_key)
            
            thetabest[i] = best
            parambest['alpha'] = best
            
        elif i == 3:
            best, output, param_true, propert = stat(i, sampler, 'b', test_params, propert)
            
            plots.stat(colours[i], output, param_true, 'beta', lnprob, zpicks, 
                       mag, sigma, nsteps, nwalkers, save_path, firstderivs_key)
            
            thetabest[i] = best
            parambest['beta'] = best
                    
        elif i == 4:
            best, output, param_true, propert = stat(i, sampler, 'g', test_params, propert)
            
            plots.stat(colours[i], output, param_true, 'gamma', lnprob, zpicks, 
                       mag, sigma, nsteps, nwalkers, save_path, firstderivs_key)
            
            thetabest[i] = best
            parambest['gamma'] = best
            
        elif i == 5:
            best, output, param_true, propert = stat(i, sampler, 'z', test_params, propert)
            
            plots.stat(colours[i], output, param_true, 'zeta', lnprob, zpicks, 
                       mag, sigma, nsteps, nwalkers, save_path, firstderivs_key)
            
            thetabest[i] = best
            parambest['zeta'] = best            

    # Checking if best found parameters are within prior.
    lp = ln.lnprior(thetabest, firstderivs_key)
    if not np.isfinite(lp):
        print('')
        print('best emcee parameters outside of prior (magbest calculation)')
        print('')

    # Plot of data mag and redshifts, overlayed with
    # mag simulated using emcee best parameters and data redshifts.
    magbest = datasim.magn(parambest, data_dict, firstderivs_key)
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
    plt.show(block=False)

    # Corner plot (walkers' walk + histogram).
    import corner
#    samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
    samples = sampler.chain[:, :, :].reshape((-1, ndim))
    corner.corner(samples, labels=["$m$", "$M$", "$alpha$", "$beta$", "$g$", "$z$"], 
                        truths=true)
    
    # Results getting printed:
    if bi == 0: 
        print('@@@@@@@@@@@@@@@@@')
        print('best index =',str(bi))
        print('@@@@@@@@@@@@@@@@@')
    print('best parameters =',str(parambest))
    print('m.a.f.:', np.mean(sampler.acceptance_fraction))
    print('nsteps:', str(nsteps))
    print('sigma:', str(sigma))
    print('npoints:', str(len(zpicks)))
    print('model:', firstderivs_key)
    
    tools.timer('burnin', timeb0, timeb1)
    tools.timer('sampler', times0, times1)
    
    return propert, sampler