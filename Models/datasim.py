3#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:47:44 2018

@author: BallBlueMeercat
"""

from math import log10
import numpy as np
import zodesolve


# Empirical parameters.
M = -19                     # Absolute brightness of supernovae.

def mag(params, zpicks):
    """
    Takes in:
            gamma = interaction constant;
            m = e_m(t)/ec(t0) at t=t0;
            de = e_de(t)/ec(t0) at t=t0;
            zpicks = list of z to match the interpolated dlmpc to;
    Returns:
        mag = list of n apparent magnitudes mag corresponding to given redshits.
    """
#    print('@@@ zmsim has been called')
        
    if not sorted(zpicks) == zpicks:
        zpicks.sort()
        print('sorted to accending in zmsim')    
    
    if not isinstance(zpicks, (list,)):
        zpicks = zpicks.tolist()
        print('converted to list in zmsim')
    
    t, dlpc, dl, a, ombar_m, ombar_de, ombar_de0 = zodesolve.zodesolve(params, zpicks)
    
    # Calculating apparent magnitudes of supernovae at the simulated
    # luminosity distances using the distance modulus formula.

    mag = []   
    for i in range(len(dlpc)):
        if dlpc[i] == 0:
            magnitude = M
        else:
            # magnitude from the distance modulus formula
            magnitude = 5 * log10(dlpc[i]/10) + M
        mag.append(magnitude)
    
#    import plots
#    de = ombar_de0
#    plots.modelcheck(t, mag, zpicks, dlpc, dl, gamma, m, de, a, ombar_m, ombar_de)
        
    return mag

def gnoise(mag, mu, sigma):
    """
   Returns:
       mag = mag, each point offset by unique Gaussian noise;
       noise = Gaussian noise.
    """
#    print('               -gnoise has been called')
    n = len(mag)
    noise = np.random.normal(mu,sigma,n)
    mag = mag + noise
    
#    import matplotlib.pyplot as pl
#    from pylab import figure
#    figure()
#    pl.title('Noise distribution')
#    pl.hist(noise, 100)
#    pl.show()
    return mag










