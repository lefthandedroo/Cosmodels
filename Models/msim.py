3#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:47:44 2018

@author: BallBlueMeercat
"""

from math import log10
import numpy as np
import odesolve

import lnprior

# Empirical parameters.
M = -19                     # Absolute brightness of supernovae.

def msim(gamma, m, de, zpicks):
    """
    Takes in:
            gamma = interaction constant;
            m = e_m(t)/ec(t0) at t=t0;
            de = e_de(t)/ec(t0) at t=t0;
            zpicks = list of z to match the interpolated dlmpc to;
    Returns:
        mag = list of n apparent magnitudes mag from corresponding redshits.
    """
#    print('@@@ msim has been called')
    
    theta = gamma, m, de
    lp = lnprior.lnprior(theta)
    if not np.isfinite(lp):
        print('msim got bad theta: ', theta)
        
    z, dlpc = odesolve.odesolve(gamma, m, de, zpicks)
    dlpcinterp = np.interp(zpicks, z, dlpc)
    
#    print('dlpcinterp is:')
#    print(dlpcinterp)

    # Calculating apparent magnitudes of supernovae at the simulated
    # luminosity distances using the distance modulus formula.
    mag = []
    for i in range(len(dlpcinterp)):
        mdistmod = 5 * log10(dlpcinterp[i]/10) + M
        mag.append(mdistmod)

    return mag