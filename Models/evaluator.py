#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 21:45:40 2018

@author: BallBlueMeercat
"""
from pylab import figure, xlabel, ylabel, title
import matplotlib.pyplot as pl
import time
import os.path

import paramfinder

# Model parameteres:  
m_true = 0.3           # (= e_m(t)/e_crit(t0) at t=t0).
de_true = 1 - m_true   # (de = e_de(t)/e_crit(t0) at t=t0).
gamma_true = 0.0       # Interaction term, rate at which DE decays into matter.

# Number of datapoints to be simulated and number of emcee steps.
npoints, nsteps = 10000, 10000

# Statistical parameteres:
mu = 0          # mean
sigma = 0.085     # standard deviation

def repeatrun():
    # Folder for saving output.
    directory = 'run'+str(int(time.time()))
    print('run directory:',directory)
    if not os.path.exists(directory):
        os.makedirs(directory)
        
    i = 0
    while i < 1:
        print('_____________________ run number',i)
        propert, sampler = paramfinder.paramfinder(npoints, nsteps, sigma, mu, m_true, directory)
        i += 1
    
    # Saving sampler to directory.
    import save
    save.output(directory, 'sampler', sampler)
    
    return sampler

#sampler = repeatrun()


def errorvsdatasize():
    # Folder for saving output.
    directory = str(int(time.time()))
    print('output directory:',directory)
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    
    sigma = 0.005
    run = 0
    
    sigma_l = []
    npoints_l = []
    sd_l = []
    mean_l = []
    vc_l = []
    sampler_l = []
    
    while sigma < 0.2:        # 0.2

        npoints = 1000  
        while npoints < 40000:#6000:   #40000
            print('_____________________ run number',run)
            propert, sampler = paramfinder.paramfinder(npoints, nsteps, sigma, mu, m_true, directory)
            sd, mean = propert
            vc = sd/mean * 100
            vc_l.append(vc)
            sd_l.append(sd)
            mean_l.append(mean)
            sigma_l.append(sigma)
            npoints_l.append(npoints)
            sampler_l.append(sampler)
            
            npoints += 2000
            run += 1
        
        sigma += 0.005
        
        # Saving plots to run directory.
        save_path = '/Users/usyd/Documents/Study/MPhil/Geraint/Cosmodels/Models/'+directory
        
    figure()
    xlabel('size of dataset')
    ylabel('standard deviation of marginalised m distribution')
    title('sd of m vs size of dataset, sd of noise = %s'%(sigma))
    pl.scatter(npoints_l, sd_l, c='m')        
    stamp = str(int(time.time()))
    filename = str(stamp)+'_sd_of_m_.png'
    filename = os.path.join(save_path, filename)
    pl.savefig(filename)
    pl.show()
    
    figure()
    xlabel('size of dataset')
    ylabel('mean of marginalised m distribution')
    title('mean of m vs size of dataset, sd of noise = %s'%(sigma))
    pl.scatter(npoints_l, mean_l, c='c')        
    stamp = str(int(time.time()))
    filename = str(stamp)+'_mean_of_m_.png'
    filename = os.path.join(save_path, filename)
    pl.savefig(filename)
    pl.show()
    
    figure()
    xlabel('size of dataset')
    ylabel('variance coefficient in %')
    title('sd/mean of m vs size of dataset, sd of noise = %s'%(sigma))
    pl.scatter(npoints_l, vc_l, c='coral')        
    stamp = str(int(time.time()))
    filename = str(stamp)+'_cv_of_m_.png'
    filename = os.path.join(save_path, filename)
    pl.savefig(filename)
    pl.show()

        
    # Saving results to directory.
    import results
    results.save(directory, 'vc', vc_l)
    results.save(directory, 'sd', sd_l)
    results.save(directory, 'mean', mean_l)
    results.save(directory, 'sigma', sigma_l)
    results.save(directory, 'npoints', npoints_l)
    results.save(directory, 'sampler', sampler_l)
    
    return vc_l, sd_l, mean_l, sigma_l, npoints_l, sampler_l,

vc, sd, mean, sigma, npoints, sampler = errorvsdatasize()