#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 16:02:10 2018

@author: BallBlueMeercat
"""

from pylab import figure, scatter, xlabel, ylabel, title, show, savefig
from pylab import legend, errorbar
from emcee import EnsembleSampler
import numpy as np
import time
import os.path

import datasim
from tools import timer
import ln
from plots import stat

def stats(params, zpicks, mag, sigma, nsteps, 
          save_path, firstderivs_key, plot_key):
    """
    Takes in:
            m_true = e_m(t)/ec(t0) at t=t0;
            zpicks = list of z to match the interpolated dlmpc to;
            mag = list of n apparent magnitudes mag for zpicks redshits;
            sigma = standard deviation used to generate Gaussian noise.
    Returns:
    """
#    print('-stats has been called')
    
    # emcee parameters:
    ndim = len(params)
    nwalkers = int(ndim * 2)

    # Initializing walkers in a Gaussian ball around the max likelihood.
    # Number in front of the np.random.rand(ndim) is 'initial footprint'.
    poslist = list(params.values())
    pos = []
    for i in poslist:
        pos.append(i / 2)
    startpos = np.array(pos)
    pos = [startpos + 0.01*np.random.randn(ndim) for i in range(nwalkers)]
    
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
                                    args=(zpicks, mag, sigma, firstderivs_key, ndim))
    
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
    for i in range(ndim):
        if i == 0:                       
            mbest = sampler.flatchain[bi,i]
            thetabest[i] = mbest
            parambest['m'] = mbest
            # Input m = e_m(z)/ec(z=0).
            m_true = params.get('m', 0)
            true.append(m_true)
            # Output m.
            m = sampler.flatchain[:,i]
            # Standard deviation and mean of the m distribution.
            m_sd = np.std(m)
            m_mean = np.mean(m)
            propert['m_sd'] = m_sd
            propert['m_mean'] = m_mean
            propert['m'] = mbest
            stat('coral', m, m_true, 'Matter', lnprob, zpicks, 
                 mag, sigma, nsteps, nwalkers, save_path)
            
        elif i == 1:
            gammabest = sampler.flatchain[bi,i]
            thetabest[i] = gammabest
            parambest['gamma'] = gammabest
            # Input interaction term.
            g_true = params.get('gamma',0)
            true.append(g_true)
            # Output gamma.
            gamma = sampler.flatchain[:,i]
            # Standard deviation and mean of the gamme distribution.
            gamma_sd = np.std(gamma)
            gamma_mean = np.mean(gamma)
            propert['gamma_sd'] = gamma_sd
            propert['gamma_mean'] = gamma_mean
            propert['gamma'] = gammabest
            stat('aquamarine', gamma, g_true, 'Gamma', lnprob, zpicks, 
                 mag, sigma, nsteps, nwalkers, save_path)
                    
        elif i == 2:
            debest = sampler.flatchain[bi,i]
            thetabest[i] = debest
            parambest['de'] = debest
            H0 = 1
            rho_c0 = H0**2 # critical density
            # Input de = e_de(z)/ec(z=0).
            de_true = params.get('de', rho_c0/rho_c0 - m_true)
            true.append(de_true)
            # Output de.
            de = sampler.flatchain[:,i]
            # Standard deviation and mean of the de distribution
            de_sd = np.std(de)
            de_mean = np.mean(de)
            propert['de_sd'] = de_sd
            propert['de_mean'] = de_mean
            propert['de'] = debest
            stat('orchid', de, de_true, 'DE', lnprob, zpicks, 
                 mag, sigma, nsteps, nwalkers, save_path)
            
    # Checking if best found parameters are within prior.
    lp = ln.lnprior(thetabest, firstderivs_key)
    if not np.isfinite(lp):
        print('')
        print('best emcee parameters outside of prior (magbest calculation)')
        print('')

    # Plot of data mag and redshifts, overlayed with
    # mag simulated using emcee best parameters and data redshifts.
    magbest = datasim.magn(parambest, zpicks, firstderivs_key, plot_key)
    figure()
    title('Evolution of magnitude with redshift \n nsteps: '
          +str(nsteps)+', noise: '+str(sigma)+', npoints: '+str(len(zpicks)))
    data = errorbar(zpicks, mag, yerr=sigma, fmt='.', alpha=0.3)
    best_fit = scatter(zpicks, magbest, lw='1', c='xkcd:tomato')
    ylabel('magnitude')
    xlabel('z')
    legend([data, best_fit], ['LCDM', firstderivs_key])
    stamp = str(int(time.time()))
    filename = str(stamp)+'____magz__nsteps_'+str(nsteps)+'_nwalkers_' \
    +str(nwalkers)+'_noise_'+str(sigma)+'_numpoints_'+str(len(zpicks))+'.png'
    filename = os.path.join(save_path, filename)
    savefig(filename)
    show()

#    # Corner plot (walkers' walk + histogram).
#    import corner
##    samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
#    samples = sampler.chain[:, :, :].reshape((-1, ndim))
#    corner.corner(samples, labels=["$m$"], 
#                        truths=true)
#    show()

    # Results getting printed:
    if bi == 0: 
        print('best index =',str(bi))
    print('best parameters =',str(parambest.values()))
    print('m.a.f.:', np.mean(sampler.acceptance_fraction))
    print('nsteps:', str(nsteps))
    print('sigma:', str(sigma))
    print('npoints:', str(len(zpicks)))
    
    timer('burnin', timeb0, timeb1)
    timer('sampler', times0, times1)
    
    return propert, sampler