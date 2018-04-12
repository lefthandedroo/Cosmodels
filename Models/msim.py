3#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:47:44 2018

@author: BallBlueMeercat
"""

from math import log10
import numpy as np
import odesolve

#import lnprior

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
    
#    theta = gamma, m, de
#    lp = lnprior.lnprior(theta)
#    if not np.isfinite(lp):
#        print('msim got bad theta: ', theta)
        
    z, dlpc, dl, gamma, e_dash0m, e_dash0de, t, a, a_dot, t_cut, a_cut, a_dotcut, e_dashm, e_dashde = odesolve.odesolve(gamma, m, de, zpicks)
    dlpcinterp = np.interp(zpicks, z, dlpc)
    
#    print('dlpcinterp is:')
#    print(dlpcinterp)

    # Calculating apparent magnitudes of supernovae at the simulated
    # luminosity distances using the distance modulus formula.
    mag = []
    for i in range(len(dlpcinterp)):
        mdistmod = 5 * log10(dlpcinterp[i]/10) + M
        mag.append(mdistmod)
    
    
    print('after msim mdistmod calculation')
    print('len dlpc is: ',len(dlpc))
    print('len dl is: ',len(dl))
    print('len a is: ',len(a))
    print('len e_dashm is: ',len(e_dashm))
    print('len e_dashde is: ',len(e_dashde))
    print('len mag is: ',len(mag))
    print('len zpicks is: ',len(zpicks))
    
#    import plots
#    plots.plots(mag, zpicks, z, dlpc, dl, gamma, e_dash0m, e_dash0de, t, a, a_dot, t_cut, a_cut, a_dotcut, e_dashm, e_dashde)
        
    theta = t, mag, dlpc, dl, a, e_dashm, e_dashde
    return theta #mag