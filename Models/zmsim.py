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
#    print('@@@ zmsim has been called')
    
    if abs(zpicks[0]) > 0:
        zpicks = [0.0] + zpicks
        print('0 added at the front')
    
    if not isinstance(zpicks, (list,)):
        zpicks = zpicks.tolist()
        print('converted to list')
   
    if not sorted(zpicks) == zpicks:
        zpicks.sort()
        print('sorted to accending')
    
#    print('len zpicks in zmsim', len(zpicks))
#    print('zpicks inside zmsim:', zpicks)
#    print(type(zpicks))
    
    dlpc, dl, gamma, e_dash0m, e_dash0de, a, e_dashm, e_dashde = zodesolve.zodesolve(gamma, m, de, zpicks)
    
    # Calculating apparent magnitudes of supernovae at the simulated
    # luminosity distances using the distance modulus formula.

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
#    print('len zpicks is: ',len(zpicks))   

#    print('dlpc is')
#    print(dlpc)
    
    mag = []
        
    for i in range(len(dlpc)):
        if dlpc[i] == 0:
            i += 1
        mdistmod = 5 * log10(dlpc[i]/10) + M
        mag.append(mdistmod)
    
#    print('after mdistmod calculation')
#    print('len dlpc is: ',len(dlpc))
#    print('len dl is: ',len(dl))
#    print('len a is: ',len(a))
#    print('len e_dashm is: ',len(e_dashm))
#    print('len e_dashde is: ',len(e_dashde))
#    print('len mag is: ',len(mag))
#    print('len zpicks is: ',len(zpicks))
    
    import zplots
    zplots.zplots(mag, zpicks, dlpc, dl, gamma, e_dash0m, e_dash0de, a, e_dashm, e_dashde)

    return #mag