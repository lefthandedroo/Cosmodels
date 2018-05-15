#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 16:02:10 2018

@author: BallBlueMeercat
"""

from pylab import figure, scatter, xlabel, ylabel, title, plot, show
import matplotlib.pyplot as pl
import emcee
import numpy as np
import time
import os.path

import tools
import ln

# emcee parameters:
ndim, nwalkers = 1, 2

def stats(m_true, zpicks, mag, sigma, nsteps, directory):
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

    # Are walkers starting outside of prior?
    i = 0
    while i < nwalkers:
        posrow = pos[i]
        theta = posrow[0]
        lp = ln.lnprior(theta)
        if not np.isfinite(lp):
            print('~~~~~~~theta %s (outside of prior) = %s ~~~~~~~'%(i, theta))
        i += 1
        
    # Sampler setup
    timee0 = time.time()    # starting sampler timer
    print('_____ sampler start')
    sampler = emcee.EnsembleSampler(nwalkers, ndim, ln.lnprob, args=(zpicks, mag, sigma))
    
    # burnin
    burnin = int(nsteps/5)     # steps to discard
    pos, prob, state = sampler.run_mcmc(pos, burnin)
    sampler.reset()
    
    # Starting sampler after burnin.
    sampler.run_mcmc(pos, nsteps)
    print('_____ sampler end')
    timee1=time.time()      # stopping sampler timer
    
    # Best parameters found by emcee.
    bi = np.argmax(sampler.flatlnprobability) # index with highest posterior prob                                       
    mbest = sampler.flatchain[bi,0]
    
    theta = mbest
    lp = ln.lnprior(theta)
    if not np.isfinite(lp):
        print('')
        print('best emcee parameters are outside of prior (magbest calcualation)')
        print('')
    
    # Saving plots to run directory.
    save_path = '/Users/usyd/Documents/Study/MPhil/Geraint/Cosmodels/Models/'+directory
    
    # Plot of mag simulated using "true" parameters, overlayed with
    # mag simulated using emcee best parameters.
    import datasim
    magbest = datasim.mag(mbest, zpicks)
    figure()
    title('Evolution of magnitude with redshift \n nsteps: '+str(nsteps)+', noise: '
          +str(sigma)+', npoints: '+str(len(zpicks)))
    data = pl.errorbar(zpicks, mag, yerr=sigma, fmt='.', alpha=0.3)
    best_fit = scatter(zpicks, magbest, lw='1', c='r')
    ylabel('magnitude')
    xlabel('z')
    pl.legend([data, best_fit], ['mag simulated with true parameters', 'mag recreated with emcee parameters'])
    stamp = str(int(time.time()))
    filename = str(stamp)+'__magVz_nsteps_'+str(nsteps)+'_nwalkers_'+str(nwalkers)+'_noise'+str(sigma)+'numpoints_'+str(len(zpicks))+'.png'
    filename = os.path.join(save_path, filename)
    pl.savefig(filename)
    show()
    show()

    # Marginalised distribution histograms.
    figure()
    xlabel(r'$\Omega_m(z=0)$')
    title('Marginalised distribution of m \n nsteps: '+str(nsteps)+', noise: '
          +str(sigma)+', npoints: '+str(len(zpicks)))
    pl.hist(sampler.flatchain[:,0], 50)
    stamp = str(int(time.time()))
    filename = str(stamp)+'__mhist_nsteps_'+str(nsteps)+'_nwalkers_'+str(nwalkers)+'_noise'+str(sigma)+'numpoints_'+str(len(zpicks))+'.png'
    filename = os.path.join(save_path, filename)
    pl.savefig(filename)
    show()
    
    # Walker steps.
    m = sampler.flatchain[:,0]
    slnprob = sampler.flatlnprobability
    
    figure()
    title('slnprob for m \n nsteps: '+str(nsteps)+', noise: '
          +str(sigma)+', npoints: '+str(len(zpicks)))
    plot(m, slnprob, '.', color='red')
    stamp = str(int(time.time()))
    filename = str(stamp)+'__steps_nsteps_'+str(nsteps)+'_nwalkers_'+str(nwalkers)+'_noise'+str(sigma)+'numpoints_'+str(len(zpicks))+'.png'
    filename = os.path.join(save_path, filename)
    pl.savefig(filename)
    show()
    
    # Chains.    
    figure()
    xlabel('step number')
    ylabel(r'$\Omega_m(z=0)$')
    title('flatChains with m_true in red \n nsteps: '+str(nsteps)+', noise: '
          +str(sigma)+', npoints: '+str(len(zpicks)))
    plot(sampler.flatchain[:,0].T, '-', color='k', alpha=0.3)
    pl.axhline(m_true, color='red')
    stamp = str(int(time.time()))
    filename = str(stamp)+'__chain_nsteps_'+str(nsteps)+'_nwalkers_'+str(nwalkers)+'_noise'+str(sigma)+'numpoints_'+str(len(zpicks))+'.png'
    filename = os.path.join(save_path, filename)
    pl.savefig(filename)
    show()

#    # Corner plot (walkers' walk + histogram).
#    print('_____ stats corner plot')
#    import corner
#    samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
#    corner.corner(samples, labels=["$m$"], 
#                        truths=[m_true])
#    pl.show()

    # Results getting printed:
#    if bi ==0:
#        print('best index =',str(bi))
    print('best index =',str(bi))
    print('mbest =',str(mbest))
    print('m.a.f.:', np.mean(sampler.acceptance_fraction))
    print('nsteps:', str(nsteps))
    print('sigma:', str(sigma))
    print('npoints:', str(len(zpicks)))
    
    # Standard deviation and mean of the m distribution
    mall = sampler.flatchain[:,0]
    sd = np.std(mall)
    mean = np.mean(mall)
    
    # Property being investigated
    propert = sd, mean
        
    tools.timer('sampler', timee0, timee1)
    
    return propert, sampler