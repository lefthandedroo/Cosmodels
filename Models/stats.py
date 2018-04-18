#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 16:02:10 2018

@author: BallBlueMeercat
"""

from pylab import figure, scatter
import corner
import emcee
import matplotlib.pyplot as pl
import numpy as np
import time

import zmsim
import lnprob
import lnprior

# emcee parameters:
ndim, nwalkers = 3, 8
nsteps = 1000
burnin = 300

def stats(gamma_true, m_true, de_true, zpicks, mag, sigma):
    """
    Takes in:
            gamma_true = interaction constant;
            m_true = e_m(t)/ec(t0) at t=t0;
            de_true = e_de(t)/ec(t0) at t=t0;
            zpicks = list of z to match the interpolated dlmpc to;
            mag = list of n apparent magnitudes mag for zpicks redshits;
            sigma = standard deviation used to generate Gaussian noise.
    Returns:
    """
#    print('-stats has been called')
        
    gamma_ml = 0.0
    m_ml = 0.3
    de_ml = 0.7
    
    startpos = np.array([gamma_ml,m_ml,de_ml])
            
    # Initializing walkers in a Gaussian ball around the max likelihood.
    # Number in front of the np.random.rand(ndim) is 'initial footprint'.
    pos = [startpos + 0.1*np.random.randn(ndim) for i in range(nwalkers)]    
#    print('pos = ',pos)    
    
    # Sampler setup
    timee0 = time.time()    # starting emcee timer
    print('_____ sampler start')
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob.lnprob, args=(zpicks, mag, sigma))
    sampler.run_mcmc(pos, nsteps)
    print('_____ sampler end')
    timee1=time.time()      # stopping emcee timer
    
    print('_____ stats corner plot')
    
    # Corner plot (walkers' walk + histogram).
    samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
    fig = corner.corner(samples, labels=["$\gamma$", "$m$", "$de$"], 
                        truths=[gamma_true, m_true, de_true])
    pl.show()
    numpoints = len(zpicks)
    fig.savefig('zz_Day'+str(time.strftime("%j"))+'_nsteps_'+
                str(nsteps)+'_'+'nwalkers_'+str(nwalkers)+'numpoints_'+
                str(numpoints)+str(time.strftime("%c"))+'.png')
    
#    # Marginalised distribution (histogram) plot.
#    figure()
#    pl.title('Marginalised distribution')
#    pl.hist(sampler.flatchain[:,0], 100)
#    pl.show()
#    
#    # Chains.    
#    figure()
#    pl.title('flatChains with gamma_true in blue')
#    pl.plot(sampler.flatchain[:,0].T, '-', color='k', alpha=0.3)
#    pl.axhline(gamma_true, color='blue')
#    pl.show
#
#    figure()
#    pl.title('flatChains with m_true in red')
#    pl.plot(sampler.flatchain[:,1].T, '-', color='k', alpha=0.3)
#    pl.axhline(m_true, color='red')
#    pl.show
#    
#    figure()
#    pl.title('flatChains with de_true in green')
#    pl.plot(sampler.flatchain[:,2].T, '-', color='k', alpha=0.3)
#    pl.axhline(de_true, color='green')
#    pl.show
    

#    print('_____ magbest calculation')
    # Simulating magnitude using best parameters found by emcee.
    bi = np.argmax(sampler.flatlnprobability)   # index with highest post prob                                       
    gammabest = sampler.flatchain[bi,0]         # parameters with the highest 
    mbest = sampler.flatchain[bi,1]             # posterior probability
    debest = sampler.flatchain[bi,2]
    

    # Results getting printed:
    print('best index is =',str(bi))
    print('gammabest is =',str(gammabest))
    print('mbest is =',str(mbest))
    print('debest is =',str(debest))
  
    # Mean acceptance fraction. In general, acceptance fraction has an entry 
    # for each walker so, in this case, it is a nwalkers-dimensional vector.
    print('Mean acceptance fraction:', np.mean(sampler.acceptance_fraction))
    print('Number of steps:', str(nsteps))
    print('Number of walkers:', str(nwalkers))
    
    theta = gammabest, mbest, debest
    lp = lnprior.lnprior(theta)
    if not np.isfinite(lp):
        print('')
        print('parameters are outside of prior when they get to magbest')
        print('')
    # Plot of mag simulted using best emcee parameters.
    magbest = zmsim.zmsim(gammabest, mbest, debest, zpicks)
    
    
    # Plot of magnitudes simulated using "true" parameters, overlayed with
    # magnitudes simulated using emcee best parameters.
    figure()
    pl.title('True parameters mag and best emcee parameters mag')
    pl.errorbar(zpicks, mag, yerr=sigma, fmt='o', alpha=0.3)
    best_fit = scatter(zpicks, magbest, lw='3', c='r')
    pl.legend([best_fit], ['Mag simulated with best emcee parameters'])
    pl.show()
        
    slnprob = sampler.flatlnprobability
    gamma = sampler.flatchain[:,0]
    m = sampler.flatchain[:,1]
    de = sampler.flatchain[:,2]
        
#    import check
#    check.thetainf(gamma, m, de, slnprob)
    
    import rslt
    rslt.rslt(gamma, m, de, slnprob)
    
    import timer
    timer.timer('sampler', timee0, timee1)
    
    return gamma, m, de, slnprob, pos