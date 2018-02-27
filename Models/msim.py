#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:47:44 2018

@author: BallBlueMeercat
"""

from math import log10
import numpy as np
import odesolve

# Empirical parameters.
M = -19                     # Absolute brightness of supernovae.

def msim(gamma, m, de, n, zpicks):
    """
    Takes in:
            gamma = interaction constant;
            m = e_m(t)/ec(t0) at t=t0;
            de = e_de(t)/ec(t0) at t=t0;
            n = dimensionless number of data points to be generated;
            zpicks = list of z to match the interpolated dlmpc to.
    Returns:
        mag = list of n apparent magnitudes mag from corresponding redshits.
    """
#    print('@@@ msim has been called')
    z, dlpc = odesolve.odesolve(gamma, m, de)
    dlpcinterp = np.interp(zpicks, z, dlpc)

    # Calculating apparent magnitudes of supernovae at the simulated
    # luminosity distances using the distance modulus formula.
    mag = []
    for i in range(len(dlpcinterp)):
        mdistmod = 5 * log10(dlpcinterp[i]/10) + M
        mag.append(mdistmod) 
    return mag
