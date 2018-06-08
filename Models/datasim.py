3#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:47:44 2018

@author: BallBlueMeercat
"""
import random
from math import log10
import numpy as np
import zodesolve
import tools

def redshift_picks(zmin, zmax, n):
    """
    Takes in:
        zmin = integer lowest redshift;
        zmax = integer highest redshift;
        n = integer number of redshifts to be generated.
    Returns:
        zpicks = list of randomly selected redshifts between zmin and zmax.
    """
#    print('-zpicks has been called')
    zinterval = (zmax - zmin) / (n*2)
    z_opts = tools.flist(zmin, zmax, zinterval)
    zpicks = random.sample(z_opts, n)
    zpicks = sorted(zpicks)
    return zpicks


def magn(params, zpicks, firstderivs_key):
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
    
    # Absolute brightness of supernovae.
    M = -19
    
    if not sorted(zpicks) == zpicks:
        zpicks.sort()
        print('sorted to accending in zmsim')    
    
    if not isinstance(zpicks, (list,)):
        zpicks = zpicks.tolist()
        print('converted to list in zmsim')
    
    dlpc, plot_var = zodesolve.zodesolve(params, zpicks, firstderivs_key)
    
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
    
#    # Checking evolution of the model.
#    import plots
#    plots.modelcheck(mag, zpicks, plot_var, firstderivs_key)
        
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


def data(mu, sigma, npoints, params, firstderivs_key):
    
    zpicks = redshift_picks(0.005, 2, npoints)
    
    model = magn(params, zpicks, firstderivs_key)
    model = np.asarray(model)
    mag = gnoise(model, mu, sigma)
    
    return mag, zpicks







