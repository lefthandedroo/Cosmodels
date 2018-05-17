#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:50:29 2017

@author: BallBlueMeercat
"""

# Model parameteres:  
m_true = 0.3           # (= e_m(t)/e_crit(t0) at t=t0).
de_true = 1 - m_true   # (de = e_de(t)/e_crit(t0) at t=t0).
gamma_true = 0.0       # Interaction term, rate at which DE decays into matter.

params = {'m_true':m_true, 'gamma_true':gamma_true}

ombar_m0 = params.get('m_true', 'not found')
rho_c0 = 1
de = rho_c0/rho_c0 - ombar_m0
print('ombar_de0 should be',de)
ombar_de0 = params.get('de_true', rho_c0/rho_c0 - ombar_m0)
print('ombar_m0',ombar_m0)
print('de',ombar_de0)



from pylab import figure, plot, xlabel, ylabel, title, show
import matplotlib.pyplot as pl

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