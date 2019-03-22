#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 16:02:10 2018

@author: BallBlueMeercat
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from emcee import EnsembleSampler
import numpy as np
import time
import os.path

import datasim
import tools
import ln
import plots

def stats(names, values, data_dict, sigma, nsteps,
          save_path, model_key, xwalkers=10, pool=None, plot=False, filename='print.txt'):
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
        model_key = string, name of model to use for model mag.
    Returns:
        propert, sampler
    """

    zpicks = data_dict.get('zpicks',0)
    mag = data_dict.get('mag',0)

    # emcee parameters:
    ndim = len(values)
    nwalkers = int(ndim * xwalkers)
    print('stats xwalkers', xwalkers)

    # Initializing walkers.
    pos = [values + 0.001*np.random.randn(ndim) for i in range(nwalkers)]

    # Are walkers starting outside of prior?
    for i in range(nwalkers):
        theta_start = pos[i]
        lp = ln.lnprior(theta_start, model_key)
        if not np.isfinite(lp):
            print('~~~~~~~pos[%s] (outside of prior) = %s ~~~~~~~'
                  %(i, theta_start))

    f=open(filename, "a+")
    f.write('inside stats pre EnsembleSampler'+'\n')

    # Sampler setup.
    times0 = time.time()    # starting sampler timer
    sampler = EnsembleSampler(nwalkers, ndim, ln.lnprob, pool=pool,
                              args=(data_dict, sigma, model_key, names))

    f.write('inside stats starting burnin'+'\n')

    # Burnin.
    burnin = int(nsteps/4)  # steps to discard
    print('_____ burnin start')
    timeb0 = time.time()    # starting burnin timer
    pos, prob, state = sampler.run_mcmc(pos, burnin)
    timeb1=time.time()      # stopping burnin timer
    print('_____ burnin end')
    sampler.reset()

    f.write('inside stats starting sampler'+'\n')

    # Starting sampler after burnin.
    print('_____ sampler start')
#    sampler.run_mcmc(pos, nsteps)
    for i, value in enumerate(sampler.sample(pos, iterations=nsteps)):
        if (i+1) % 100 == 0:
            print("{0:5.1%}".format(float(i) / nsteps))
            f.write("{0:5.1%}".format(float(i) / nsteps) +'\n')
    print('_____ sampler end')
    times1=time.time()      # stopping sampler timer

    f.write('inside stats sampler finished'+'\n')
    f.close()

    # Walker steps is lnprob = sampler.flatlnprobability
    # Index of parameter values with highest post prob found by emcee.
    bi = np.argmax(sampler.flatlnprobability)

    # Extracting results:
    thetabest = sampler.flatchain[bi,:]
    # Checking if best found parameters are within prior.
    lp = ln.lnprior(thetabest, model_key)
    if not np.isfinite(lp):
        print('~~~~~~~thetabest outside of prior at magbest~~~~~~~')

    propert = {}
    propert['trace'] = sampler.chain[:, burnin:, :].reshape(-1, ndim) #trace

    colours = ['light red','berry', 'coral', 'amber', 'apple',
        'aquamarine', 'raspberry', 'green blue', 'deep blue',
        'emerald', 'blue violet', 'dark violet', 'yellow orange',
        'light red', 'berry', 'coral', 'amber', 'apple', 'aquamarine',
        'raspberry', 'green blue', 'deep blue', 'emerald', 'blue violet',
        'dark violet', 'black']

    for i in range(ndim):
        param_initial = names[i][0]
        # Fitted parameter.
        output = sampler.flatchain[:,i]
        # Standard deviation and mean of the emcee found distribution.
        propert[param_initial+'_sd'] = np.std(output)
        propert[param_initial+'_mean'] = np.mean(output)
        propert[param_initial] = sampler.flatchain[bi,i]

        if plot:
            plots.stat_emcee(colours[i], output, values[i], names[i],
                       sampler.flatlnprobability, zpicks, mag, sigma,
                       nsteps, nwalkers, save_path, model_key)

    if plot:
        # Plot of data mag and redshifts, overlayed with
        # mag simulated using emcee best parameters and data redshifts.
        magbest, dabest = datasim.magn(names, thetabest, data_dict, model_key)
        plt.figure()
        plt.title('model: '+model_key
            +'\n Evolution of magnitude with redshift \n nsteps: '
            +str(nsteps)+', noise: '+str(sigma)+', npoints: '+str(len(zpicks)))
        data = plt.errorbar(zpicks, mag, yerr=sigma, fmt='.', alpha=0.3)
        best_fit = plt.scatter(zpicks, magbest, lw='1', c='xkcd:tomato')
        plt.ylabel('magnitude')
        plt.xlabel('z')
        plt.legend([data, best_fit], ['LCDM', model_key])
        stamp = str(int(time.time()))
        filename = str(stamp)+'____magz__nsteps_'+str(nsteps)+'_nwalkers_' \
        +str(nwalkers)+'_noise_'+str(sigma)+'_numpoints_'+str(len(zpicks))+'.png'
        filename = os.path.join(save_path, filename)
        plt.savefig(filename)

        # Corner plot (walkers' walk + histogram).
        import corner
        mpl.style.use('default')
        samples = sampler.chain[:, :, :].reshape((-1, ndim))
        corner.corner(samples, labels=names, truths=values)
        stamp = str(int(time.time()))
        filename = str(stamp)+'____corn__nsteps_'+str(nsteps)+'_nwalkers_' \
        +str(nwalkers)+'_noise_'+str(sigma)+'_numpoints_'+str(len(zpicks))+'.png'
        filename = os.path.join(save_path, filename)
        plt.savefig(filename)

        plt.show(block=False)

    # Results getting printed:
    if bi == 0:
        print('@@@@@@@@@@@@@@@@@')
        print('best index =',str(bi))
        print('@@@@@@@@@@@@@@@@@')
    print('max likelihood params =',str(thetabest))
    print('m.a.f.:', np.mean(sampler.acceptance_fraction))
    print('nsteps:', str(nsteps))
    print('ndim:',str(ndim))
    print('nwalkers:', str(nwalkers))
    print('len(sampler.chain[:, :, :].reshape((-1, ndim))):', str(len(sampler.chain[:, :, :].reshape((-1, ndim)))))
    print('sigma:', str(sigma))
    print('npoints:', str(len(zpicks)))
    print('model:', model_key)

    tools.timer('burnin', timeb0, timeb1)
    tools.timer('sampler', times0, times1)

    return propert, sampler