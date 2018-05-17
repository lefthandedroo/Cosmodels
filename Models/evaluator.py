#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 21:45:40 2018

@author: BallBlueMeercat
"""
from pylab import figure, xlabel, ylabel, title, scatter, show, savefig
import time
import os.path

from results import save
import paramfinder
import tools

# Model parameteres:  
m_true = 0.3           # (= e_m(t)/e_crit(t0) at t=t0).
de_true = 1 - m_true   # (de = e_de(t)/e_crit(t0) at t=t0).
gamma_true = 0.0       # Interaction term, rate at which DE decays into matter.

params = {'m_true':m_true}#, 'de_true':de_true, 'gamma_true':gamma_true}

# Number of datapoints to be simulated and number of emcee steps.
npoints, nsteps = 10000, 10000

# Statistical parameteres:
mu = 0          # mean
sigma = 0.085     # standard deviation

def repeatrun():
    nsteps = 1000
    # Folder for saving output.
    directory = 'run'+str(int(time.time()))
    # Relative path of output folder.
    save_path = './'+directory 
    if not os.path.exists(save_path):
        os.makedirs(save_path)
   
    
    i = 0
    while i < 1:
        print('_____________________ run number',i)
        propert, sampler = paramfinder.paramfinder(
                npoints, nsteps, sigma, mu, params, save_path)
        i += 1
        sd, mean = propert
        
    figure()
    title('test of directory options')
    scatter(sd, mean, c='m')        
    stamp = str(int(time.time()))
    filename = str(stamp)+'_sd_of_m_.png'
    filename = os.path.join(save_path, filename)
    savefig(filename)
    show()
    
    # Saving sampler to directory.
    save(save_path, 'sampler', sampler)
    
    print('run directory:',directory)

    return sampler

sampler = repeatrun()


def errorvsdatasize():
    # Script timer.
    timet0 = time.time()
    
    # Folder for saving output.
    directory = str(int(time.time()))
    # Relative path of output folder.
    save_path = './results/'+directory
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    
    sigma = 0.15
    run = 0
    
    sigma_l = []
    npoints_l = []
    sd_l = []
    mean_l = []
    vc_l = []
    sampler_l = []
    
    while sigma < 0.2:

        npoints = 4500 
        while npoints < 10000:
            print('_____________________ run number',run)
            propert, sampler = paramfinder.paramfinder(
                    npoints, nsteps, sigma, mu, m_true, save_path)
            sd, mean = propert
            vc = sd/mean * 100
            vc_l.append(vc)
            sd_l.append(sd)
            mean_l.append(mean)
            sigma_l.append(sigma)
            npoints_l.append(npoints)
            sampler_l.append(sampler)
            
            npoints += 3000
            run += 1
        
        sigma += 0.03
        
        # Saving plots to run directory.
    figure()
    xlabel('size of dataset')
    ylabel('standard deviation of marginalised m distribution')
    title('sd of m vs size of dataset, sd of noise = %s'%(sigma))
    scatter(npoints_l, sd_l, c='m')        
    stamp = str(int(time.time()))
    filename = str(stamp)+'_sd_of_m_.png'
    filename = os.path.join(save_path, filename)
    savefig(filename)
    show()
    
    figure()
    xlabel('size of dataset')
    ylabel('mean of marginalised m distribution')
    title('mean of m vs size of dataset, sd of noise = %s'%(sigma))
    scatter(npoints_l, mean_l, c='c')        
    stamp = str(int(time.time()))
    filename = str(stamp)+'_mean_of_m_.png'
    filename = os.path.join(save_path, filename)
    savefig(filename)
    show()
    
    figure()
    xlabel('size of dataset')
    ylabel('variance coefficient in %')
    title('sd/mean of m vs size of dataset, sd of noise = %s'%(sigma))
    scatter(npoints_l, vc_l, c='coral')        
    stamp = str(int(time.time()))
    filename = str(stamp)+'_cv_of_m_.png'
    filename = os.path.join(save_path, filename)
    savefig(filename)
    show()

        
    # Saving results to directory.
    save(save_path, 'vc', vc_l)
    save(save_path, 'sd', sd_l)
    save(save_path, 'mean', mean_l)
    save(save_path, 'sigma', sigma_l)
    save(save_path, 'npoints', npoints_l)
    save(save_path, 'sampler', sampler_l)
    
    print('output directory:',directory)
    
    # Time taken by evaluator. 
    timet1=time.time()
    tools.timer('evaluator', timet0, timet1)
    
    return vc_l, sd_l, mean_l, sigma_l, npoints_l, sampler_l,

#vc, sd, mean, sigma, npoints, sampler = errorvsdatasize()