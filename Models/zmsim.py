3#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:47:44 2018

@author: BallBlueMeercat
"""

from math import log10
import zodesolve


# Empirical parameters.
M = -19                     # Absolute brightness of supernovae.

def zmsim(gamma, m, zpicks):
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
    
    t, dlpc, dl, a, ombar_m, ombar_de, ombar_de0 = zodesolve.zodesolve(gamma, m, zpicks)
    
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
    
#    import zplots
#    de = ombar_de0
#    zplots.zplots(t, mag, zpicks, dlpc, dl, gamma, m, de, a, ombar_m, ombar_de)
    
#    ztheta = t, mag, dlpc, dl, a, ombar_m, ombar_de
    
    return mag #ztheta












