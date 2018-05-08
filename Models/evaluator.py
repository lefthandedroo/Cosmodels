#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 21:45:40 2018

@author: BallBlueMeercat
"""
from pylab import figure, plot, xlabel, ylabel, title, show
import matplotlib.pyplot as pl


import paramfinder

# Model parameteres:  
m_true = 0.3           # (= e_m(t)/e_crit(t0) at t=t0).
de_true = 1 - m_true   # (de = e_de(t)/e_crit(t0) at t=t0).
gamma_true = 0.0       # Interaction term, rate at which DE decays into matter.

# Number of datapoints to be simulated and number of emcee steps.
npoints, nsteps = 10000, 20000

# Statistical parameteres:
mu = 0          # mean
sigma = 0.085     # standard deviation

def repeatrun():    
    i = 0
    while i < 1:
        print('_____________________ run number',i)
        propert, sampler = paramfinder.paramfinder(npoints, nsteps, sigma, mu, m_true)
        i += 1
    
    return sampler

sampler = repeatrun()


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


def errorvsdatasize():
    error = []
    meanlist = []
    cvlist = []
    
    sigma = 0.1
    run = 0
    while sigma > 1e-9:
        numpoints = []
        standev = []
        npoints = 4000
        while npoints < 25000:    #35000
            print('_____________________ run number',run)
            error.append(sigma)
            numpoints.append(npoints)
            propert, sampler = paramfinder.paramfinder(npoints, nsteps, sigma, mu, m_true)
            sd, mean = propert
            standev.append(sd)
            meanlist.append(mean)
            cv = sd/mean
            cvlist.append(cv)
            npoints += 1500
            run += 1
        figure()
        xlabel('size of dataset')
        ylabel('standard deviation')
        plot(numpoints, standev, '.')
        title('sd of m found vs size of dataset used, sd of noise added = %s'%(sigma))
        show()
        sigma -= 0.003
    
#    # Create two subplots sharing y axis
#    fig, (noise, dataset) = pl.subplots(2, sharey=True)
#    
#    noise.plot(error, standev, 'b.')
#    noise.set(title='', ylabel='sd of noise')
#    
#    dataset.plot(numpoints, standev, 'r.')
#    dataset.set(xlabel='standard deviation of parameter distribution', ylabel='size of dataset')
#    
#    pl.show()
        
#    # Create two subplots sharing x axis
#    fig, (noise, dataset) = pl.subplots(2, sharex=True)
#    
#    noise.plot(standev, error, 'b.')
#    noise.set(title='', ylabel='standard deviation of noise')
#    
#    dataset.plot(standev, numpoints, 'r.')
#    dataset.set(xlabel='sd of parameter distribution', ylabel='size of dataset')
#    
#    pl.show()
    
    return error, numpoints, sampler, standev, meanlist

#error, numpoints, sampler, standev, meanlist = errorvsdatasize()