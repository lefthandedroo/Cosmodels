#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main code
"""
import numpy as np
import time

import timer
import gnoise
import zmsim
import stats


## Number of datapoints to be simulated.
#n = 1000
#
## Model parameteres:  
#gamma_true = 0  # Interaction term, rate at which DE decays into matter.
#m_true = 0.3    # (= e_m(t)/e_crit(t0) at t=t0).
#de_true = 0.7   # (de = e_de(t)/e_crit(t0) at t=t0).
#
## Statistical parameteres:
#mu = 0          # mean
#sigma = 0.1     # standard deviation

def paramfinder(n, gamma_true, m_true, sigma, mu):
    # Script timer.
    timet0 = time.time()
    
    # Generating redshifts to simulate data.
    import zpicks
    zmin, zmax = 0.001, 2
    zpicks = zpicks.zpicks(zmin, zmax, n)
    zpicks = sorted(zpicks)
    # Inserting 0 at the front to allow initial conditions.
    if abs(zpicks[0]) > 0:
        zpicks = [0.0] + zpicks
    
#    print('zpicks')
#    print(zpicks)
    
    # Generating apparent magnitues mag at redshifts z=0 to z < zmax.
    model = zmsim.zmsim(gamma_true, m_true, zpicks)
    model = np.asarray(model)
    mag = gnoise.gnoise(model, mu, sigma)
    
    # Stats timer.
    times0 = time.time()
    
#    gammabest, mbest, debest = stats.stats(gamma_true, m_true, zpicks, mag, sigma)
    gammabest, mbest = stats.stats(gamma_true, m_true, zpicks, mag, sigma)

    
    # Time taken by stats. 
    times1=time.time()
    timer.timer('stats', times0, times1)
    
    # Time taken by the script. 
    timet1=time.time()
    timer.timer('script', timet0, timet1)
    
    return gammabest, mbest#, debest

#paramfinder(n, gamma_true, m_true, de_true, sigma, mu)