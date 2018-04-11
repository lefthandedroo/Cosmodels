7#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:38:48 2018

@author: BallBlueMeercat
"""
import numpy as np
from scipy.integrate import odeint

import zfirstderivs


# Standard cosmological parameters.
H0 = 1
tH = 1.0/H0  # Hubble time
c_over_H0 = 4167 * 10**6    # c/H0 in parsecs

def zodesolve(gamma,m,de,zpicks):
    """
    Takes in:
        gamma = interaction constant;
        m = e_m(t)/ec(t0) at t=t0;
        de = e_de(t)/ec(t0) at t=t0.
    Returns: 
        z = numpoints number of redshifts zmin<z<zmax;
        dlpc = luminosity distance in pc.
    
    """
#    print('@@ zodesolve has been called')
    
    # Initial conditions at z = 0.
    a0 = 1.0        # scale factor
#    a_dz0 = 1.0    # speed of expansion
    e_dash0m = m    # e_m(t)/ec(t0)
    e_dash0de = de  # e_de(t)/ec(t0)
    z0 = 0
    dl0 = 0
    
    # ODE solver parameters:
    abserr = 1.0e-8
    relerr = 1.0e-6
    
    # Pack up the initial conditions and eq of state parameters.
    v0 = [a0, e_dash0m, e_dash0de, z0, dl0]
    
    # Call the ODE solver. maxstep=5000000 added later to try and avoid 
    vsol = odeint(zfirstderivs.zfirstderivs, v0, zpicks, args=(gamma,), 
                  atol=abserr, rtol=relerr, mxstep=5000000)
            
    # Remove unwanted results which are too close to big bang from the plot.
    # Separate results into their own arrays:
    a = vsol[:,0]
    e_dashm = vsol[:,1]
    e_dashde = vsol[:,2]
    z = vsol[:,3]    
    dl = vsol[:,4] * (1+z)   # in units of dl*(H0/c)
    dlpc = dl * c_over_H0    # dl in parsecs (= vsol[dl] * c/H0)

        
#    print('len zpicks after vsol',len(zpicks))
    
    
    return dlpc, dl, gamma, e_dash0m, e_dash0de, a, e_dashm, e_dashde
