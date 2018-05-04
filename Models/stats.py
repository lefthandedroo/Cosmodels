#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 16:02:10 2018

@author: BallBlueMeercat
"""

from pylab import figure, scatter
import matplotlib.pyplot as pl
import emcee
import numpy as np
import time

import lnprob
import lnprior

# emcee parameters:
ndim, nwalkers = 1, 2

def stats(m_true, zpicks, mag, sigma, nsteps):
    """
    Takes in:
            m_true = e_m(t)/ec(t0) at t=t0;
            zpicks = list of z to match the interpolated dlmpc to;
            mag = list of n apparent magnitudes mag for zpicks redshits;
            sigma = standard deviation used to generate Gaussian noise.
    Returns:
    """
#    print('-stats has been called')
        
    # Initializing walkers in a Gaussian ball around the max likelihood.
    # Number in front of the np.random.rand(ndim) is 'initial footprint'.
    m_start = m_true / 2    
    startpos = np.array([m_start])
    pos = [startpos + 0.01*np.random.randn(ndim) for i in range(nwalkers)]

    i = 0
    while i < nwalkers:
        posrow = pos[i]
        m = posrow[0]
        theta = m
        lp = lnprior.lnprior(theta)
        if not np.isfinite(lp):
            print('~~~~~~~theta %s (outside of prior) = %s ~~~~~~~'%(i, theta))
        i += 1
       
    # Sampler setup
    timee0 = time.time()    # starting sampler timer
    print('_____ sampler start')
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob.lnprob, args=(zpicks, mag, sigma))
    sampler.run_mcmc(pos, nsteps)
    print('_____ sampler end')
    timee1=time.time()      # stopping sampler timer
    
    # Best parameters found by emcee.
    bi = np.argmax(sampler.flatlnprobability) # index with highest posterior prob                                       
    mbest = sampler.flatchain[bi,0]
    
    theta = mbest
    lp = lnprior.lnprior(theta)
    if not np.isfinite(lp):
        print('')
        print('parameters are outside of prior when they get to magbest')
        print('')
        
    # Plot of mag simulated using "true" parameters, overlayed with
    # mag simulated using emcee best parameters.
    import zmsim
    magbest = zmsim.zmsim(mbest, zpicks)
    figure()
    pl.title('True parameters mag and best emcee parameters mag')
    pl.errorbar(zpicks, mag, yerr=sigma, fmt='.', alpha=0.3)
    best_fit = scatter(zpicks, magbest, lw='1', c='r')
    pl.legend([best_fit], ['Mag simulated with best emcee parameters'])
    pl.show()

    # Marginalised distribution histograms.
    figure()
    pl.title('Marginalised distribution for m')
    pl.hist(sampler.flatchain[:,0], 150)
    pl.show()

    # Corner plot (walkers' walk + histogram).
    print('_____ stats corner plot')
    burnin = int(nsteps/10) # steps to discard
    import corner
    samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
    corner.corner(samples, labels=["$m$"], 
                        truths=[m_true])
    pl.show()

    # Results getting printed:
    print('best index is =',str(bi))
    print('mbest is =',str(mbest))
    print('Mean acceptance fraction:', np.mean(sampler.acceptance_fraction))
    print('nsteps:', str(nsteps))
    print('nwalkers:', str(nwalkers))
    print('npoints:', str(len(zpicks)))
        
    import timer
    timer.timer('sampler', timee0, timee1)
    
    return mbest, sampler