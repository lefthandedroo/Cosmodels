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
ndim, nwalkers = 2, 4

def stats(gamma_true, m_true, zpicks, mag, sigma, nsteps):
    """
    Takes in:
            gamma_true = interaction constant;
            m_true = e_m(t)/ec(t0) at t=t0;
#            de_true = e_de(t)/ec(t0) at t=t0;
            zpicks = list of z to match the interpolated dlmpc to;
            mag = list of n apparent magnitudes mag for zpicks redshits;
            sigma = standard deviation used to generate Gaussian noise.
    Returns:
    """
#    print('-stats has been called')
    
    burnin = int(nsteps/10)
        
    gamma_start = gamma_true / 2
    m_start = m_true / 2
    
    startpos = np.array([gamma_start,m_start])
            
    # Initializing walkers in a Gaussian ball around the max likelihood.
    # Number in front of the np.random.rand(ndim) is 'initial footprint'.
    pos = [startpos + 0.01*np.random.randn(ndim) for i in range(nwalkers)]
#    print('pos = ',pos) 

    i = 0
    while i < 4:
        posrow = pos[i]
        gamma = posrow[0]
        m = posrow[1]
        theta = gamma, m
        lp = lnprior.lnprior(theta)
        if not np.isfinite(lp):
            print('~~~~~~~theta %s (outside of prior) = %s ~~~~~~~'%(i, theta))
        i += 1
       
    
    # Sampler setup
    timee0 = time.time()    # starting emcee timer
    print('_____ sampler start')
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob.lnprob, args=(zpicks, mag, sigma))
    sampler.run_mcmc(pos, nsteps)
    print('_____ sampler end')
    timee1=time.time()      # stopping emcee timer
    
    print('_____ stats corner plot')
    
#    # Corner plot (walkers' walk + histogram).
#    import corner
#    samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
#    fig = corner.corner(samples, labels=["$\gamma$", "$m$", "$de$"], 
#                        truths=[gamma_true, m_true])
#    pl.show()
    
    # Naming of file
#    numpoints = len(zpicks)
#    fig.savefig('zz_Day'+str(time.strftime("%j"))+'_nsteps_'+
#                str(nsteps)+'_'+'nwalkers_'+str(nwalkers)+'numpoints_'+
#                str(numpoints)+str(time.strftime("%c"))+'.png')
    
    # Marginalised distribution (histogram) plot.
    figure()
    pl.title('Marginalised distribution for gamma')
    pl.hist(sampler.flatchain[:,0], 100)
    pl.show()
    
    # Marginalised distribution (histogram) plot.
    figure()
    pl.title('Marginalised distribution for m')
    pl.hist(sampler.flatchain[:,1], 100)
    pl.show()
    
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
    

    # Best parameters found by emcee.
    bi = np.argmax(sampler.flatlnprobability)   # index with highest post prob                                       
    gammabest = sampler.flatchain[bi,0]         # parameters with the highest 
    mbest = sampler.flatchain[bi,1]             # posterior probability
    debest = 1 - mbest

    # Results getting printed:
    print('best index is =',str(bi))
    print('gammabest is =',str(gammabest))
    print('mbest is =',str(mbest))
    print('1 - mbest =',str(debest))
  
    # Mean acceptance fraction. In general, acceptance fraction has an entry 
    # for each walker so, in this case, it is a nwalkers-dimensional vector.
    print('Mean acceptance fraction:', np.mean(sampler.acceptance_fraction))
    print('Number of steps:', str(nsteps))
    print('Number of walkers:', str(nwalkers))
    print('Number of datapoints used:', str(len(zpicks)))
    
    theta = gammabest, mbest#, debest
    lp = lnprior.lnprior(theta)
    if not np.isfinite(lp):
        print('')
        print('parameters are outside of prior when they get to magbest')
        print('')
    
    # Plot of magnitudes simulated using "true" parameters, overlayed with
    # magnitudes simulated using emcee best parameters.
#    import zmsim
#    magbest = zmsim.zmsim(gammabest, mbest, zpicks)

#    magerr = mag * sigma
#    plmagerr = magerr[1:]
#    figure()
#    pl.title('plmagerr distribution')
#    pl.hist(plmagerr, 100)
#    pl.show()

#    figure()
#    magerr = mag * sigma
#    pl.title('True parameters mag and best emcee parameters mag')
#    pl.errorbar(zpicks, mag, yerr=magerr, fmt='.', alpha=0.3)
#    best_fit = scatter(zpicks, magbest, lw='1', c='r')
#    pl.legend([best_fit], ['Mag simulated with best emcee parameters'])
#    pl.show()
        
#    slnprob = sampler.flatlnprobability
#    gamma = sampler.flatchain[:,0]
#    m = sampler.flatchain[:,1]
#    de = sampler.flatchain[:,2]
#    de = 1 - m
    
#    print('gamma')
#    print(len(gamma))
    
#    figure()
#    pl.title('complete ombar_m vs gamma')
#    pl.xlabel('m')
#    pl.ylabel('gamma')
#    scatter(m, gamma, lw='1')
#    pl.show()
    
#    gammacut = gamma[burnin:]
#    mcut = m[burnin:]
    
#    print('gammacut')
#    print(len(gammacut))
#    
#    figure()
#    pl.title('cut ombar_m vs gamma')
#    pl.xlabel('m')
#    pl.ylabel('gamma')
#    scatter(mcut, gammacut, lw='1')
#    pl.show()
    
#    figure()
#    pl.title('ombar_m vs ombar_de')
#    pl.xlabel('m')
#    pl.ylabel('de')
#    scatter(m, de, lw='1', c='y')
#    pl.show()
        
#    import check
#    check.thetainf(gamma, m, de, slnprob)
    
#    import rslt
#    rslt.rslt(gamma, m, de, slnprob)
    
    import timer
    timer.timer('sampler', timee0, timee1)
    
    return gammabest, mbest, debest