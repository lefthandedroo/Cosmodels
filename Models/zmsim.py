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

def zmsim(gamma, m, de, zpicks):
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
        
    z, dlpc, dl, gamma, e_dash0m, e_dash0de, t, a, a_dot, t_cut, a_cut, a_dotcut, e_dashm, e_dashde = zodesolve.zodesolve(gamma, m, de, zpicks)
    
    # Calculating apparent magnitudes of supernovae at the simulated
    # luminosity distances using the distance modulus formula.
    mag = []
        
    for i in range(len(dlpc)):
        mdistmod = 5 * log10(dlpc[i]/10) + M
        mag.append(mdistmod)
    
    zpicks = zpicks[1:]
    
#    print('len z is: ',len(z))
#    print('len dlpc is: ',len(dlpc))
#    print('len dl is: ',len(dl))
#    print('len t is: ',len(t))
#    print('len a is: ',len(a))
#    print('len a_dot is: ',len(a_dot))
#    print('len t_cut is: ',len(t_cut))
#    print('len a_cut is: ',len(a_cut))
#    print('len a_dotcut is: ',len(a_dotcut))
#    print('len e_dashm is: ',len(e_dashm))
#    print('len e_dashde is: ',len(e_dashde))
#    print('len mag is: ',len(mag))
#    print('len zpicks is: ',len(zpicks))
    
    import plots
    plots.plots(mag, zpicks, z, dlpc, dl, gamma, e_dash0m, e_dash0de, t, a, a_dot, t_cut, a_cut, a_dotcut, e_dashm, e_dashde)

    return #mag