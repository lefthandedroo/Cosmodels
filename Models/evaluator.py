#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 21:45:40 2018

@author: BallBlueMeercat
"""
from pylab import figure, plot, xlabel, ylabel, title, show
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
    
    error = []
    
    sigma = 1e-9
    run = 0
    
    while sigma < 0.006:
        numpoints = []
        standev = []
        meanlist = []
        cvlist = []
        samplerlist = []
        npoints = 4000
        while npoints < 50000:    #35000
            print('_____________________ run number',run)
            error.append(sigma)
            numpoints.append(npoints)
            propert, sampler = paramfinder.paramfinder(npoints, nsteps, sigma, mu, m_true, directory)
            sd, mean = propert
            standev.append(sd)
            meanlist.append(mean)
            cv = sd/mean * 100
            cvlist.append(cv)
            samplerlist.append(sampler)
            
            npoints += 2000
            run += 1
        
        # Saving plots to run directory.
        save_path = '/Users/usyd/Documents/Study/MPhil/Geraint/Cosmodels/Models/'+directory
        
        figure()
        xlabel('size of dataset')
        ylabel('standard deviation')
        title('sd of m vs size of dataset, sd of noise = %s'%(sigma))
        pl.scatter(numpoints, standev, c='m')        
        stamp = str(int(time.time()))
        filename = str(stamp)+'_sd_of_m_.png'
        filename = os.path.join(save_path, filename)
        pl.savefig(filename)
        pl.show()
        
        figure()
        xlabel('size of dataset')
        ylabel('mean of marginalised distribution for m')
        title('mean of m vs size of dataset, sd of noise = %s'%(sigma))
        pl.scatter(numpoints, meanlist, c='c')        
        stamp = str(int(time.time()))
        filename = str(stamp)+'_mean_of_m_.png'
        filename = os.path.join(save_path, filename)
        pl.savefig(filename)
        pl.show()
        
        figure()
        xlabel('size of dataset')
        ylabel('variance coefficient in %')
        title('sd/mean of m vs size of dataset, sd of noise added = %s'%(sigma))
        pl.scatter(numpoints, cvlist, c='coral')        
        stamp = str(int(time.time()))
        filename = str(stamp)+'_cv_of_m_.png'
        filename = os.path.join(save_path, filename)
        pl.savefig(filename)
        pl.show()
        
        sigma += 0.003
        
    # Saving results to directory.
    import save
    save.output(directory, 'error', error)
    save.output(directory, 'numpoints', numpoints)
    save.output(directory, 'sampler', sampler)
    save.output(directory, 'standev', standev)
    save.output(directory, 'meanlist', meanlist)
    return error, numpoints, sampler, standev, meanlist

error, numpoints, sampler, standev, meanlist = errorvsdatasize()


def stepevaluator():
    
    steps = []
    standev = []
    meanlist = []
    cvlist = []
    
    nsteps = 100
    run = 0
    while nsteps < 2000:
        print('_____________________ run number',run)
        steps.append(nsteps)
        propert, sampler = paramfinder.paramfinder(npoints, nsteps, sigma, mu, m_true)
        sd, mean = propert
        standev.append(sd)
        meanlist.append(mean)
        cv = sd/mean
        cvlist.append(cv)
        
        cv = sd / mean     # Coefficient of variation.
#        print('cv:',str(cv))
#        if cv < 0.008:
#            print('nsteps', nsteps)
#            break
        
        nsteps += 50
        run += 1
    
    figure()
    xlabel('emcee steps')
    ylabel('standard deviation')
    plot(steps, standev, '.')
    title('standard deviation of m found vs steps taken')
    show()
    
    return steps, sampler, standev, meanlist

#steps, sampler, standev, meanlist = stepevaluator()


def nevaluator():
    
    numpoints = []
    standev = []
    meanlist = []
    cvlist = []
    
    npoints = 100
    run = 0
    while npoints < 35000:    #35000
        print('_____________________ run number',run)
        numpoints.append(npoints)
        propert, sampler = paramfinder.paramfinder(npoints, nsteps, sigma, mu, m_true)
        sd, mean = propert
        standev.append(sd)
        meanlist.append(mean)
        cv = sd/mean
        cvlist.append(cv)
        npoints += 1000
        run += 1

    
    figure()
    xlabel('number of datapoints used')
    ylabel('standard deviation')
    plot(numpoints, standev, 'b.', label='standard deviation')
    plot(numpoints, cvlist, 'r.', label='coefficient of variation')
    pl.legend()
    title('sd and cv of m found vs dataset size')
    show()
    
    return numpoints, sampler, standev, meanlist

#numpoints, sampler, standev, meanlist = nevaluator()


def errevaluator():
    
    error = []
    standev = []
    meanlist = []
    cvlist = []
    
    sigma = 1e-9
    run = 0
    while sigma < 0.1:
        print('_____________________ run number',run)
        error.append(sigma)
        propert, sampler = paramfinder.paramfinder(npoints, nsteps, sigma, mu, m_true)
        sd, mean = propert
        standev.append(sd)
        meanlist.append(mean)
        cv = sd/mean
        cvlist.append(cv)
        sigma += 0.003
        run += 1
    
    figure()
    xlabel('standard deviation of noise')
    ylabel('standard deviation of parameter distribution')
    plot(error, standev, 'b.', label='standard deviation')
    plot(error, cvlist, 'r.', label='coefficient of variation')
    pl.legend()
    title('sd and cv of m found vs sd of noise added')
    show()
    
    return error, sampler, standev, meanlist

#error, sampler, standev, meanlist = errevaluator()