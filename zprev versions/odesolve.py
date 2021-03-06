7#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:38:48 2018

@author: BallBlueMeercat
"""
import numpy as np
from scipy.integrate import odeint

import firstderivs
#import lnprior



# Standard cosmological parameters.
H0 = 1       # Hubble parameter at t=now
tH = 1.0/H0  # Hubble time
c_over_H0 = 4167 * 10**6    # c/H0 in parsecs

def odesolve(gamma,m,de,zpicks):
    """
    Takes in:
        gamma = interaction constant;
        m = e_m(t)/ec(t0) at t=t0;
        de = e_de(t)/ec(t0) at t=t0.
    Returns: 
        z = numpoints number of redshifts zmin<z<zmax;
        dlpc = luminosity distance in pc.
    
    """
#    print('@@ odesolve has been called')
    
    # Time (in 1/H0) to integrate until.  If this time isn't long 
    # enough for 'a' to decrease to a_d then stoptime will be extended 
    # by time until a_d is reached.
    # 0.665 matter only, 0.96 standard m+de
    time = 0.8
    
    
    # Initial conditions at time = t0.
    a0 = 1.0        # scale factor
    a_dot0 = 1.0    # speed of expansion
    e_dash0m = m    # e_m(t)/ec(t0)
    e_dash0de = de  # e_de(t)/ec(t0)
    z0 = 0
    dl0 = 0
    
    # ODE solver parameters:
    abserr = 1.0e-8
    relerr = 1.0e-6
    numpoints = 100
    
    stoptime = 0 # Integrating back in time as time now is t0.
    
    z = np.array([0])
    
    stoptime -= time
        
#    theta = gamma, m, de
#    lp = lnprior.lnprior(theta)
#    if not np.isfinite(lp):
#        time += 500

    if time > 0.9:
        print('time in odesolve is: %s, gamma = %s, m = %s, de = %s'
              %(time, gamma, m, de))
    # Create time samples for the ODE solver.
    t = [stoptime * tH * float(i) / (numpoints - 1) for i in range(numpoints)]

    # Pack up the initial conditions and eq of state parameters.
    v0 = [a0, a_dot0, e_dash0m, e_dash0de, z0, dl0]
    
    # Call the ODE solver. maxstep=5000000 added later to try and avoid 
    vsol = odeint(firstderivs.firstderivs, v0, t, args=(gamma,), 
                  atol=abserr, rtol=relerr, mxstep=5000000)
            
    # Remove unwanted results which are too close to big bang from the plot.
    # Separate results into their own arrays:
    a = vsol[:,0]
    a_dot = vsol[:,1]
    e_dashm = vsol[:,2]
    e_dashde = vsol[:,3]
    z = vsol[:,4]
    dl = vsol[:,5] * (1+z)   # in units of dl*(H0/c)
    dlpc = dl * c_over_H0    # dl in parsecs (= vsol[dl] * c/H0)
    

    # Remove values after the index of first instance of z > 2.
    t_cut = np.asarray(t)
    
    a_cut = a
    a_dotcut = a_dot
    
#    # Age of the universe.
#    age = t_cut[np.argmin(t_cut)]
#    age = -round(age, 2)
    
    
    return z, dlpc, dl, gamma, e_dash0m, e_dash0de, t, a, a_dot, t_cut, a_cut, a_dotcut, e_dashm, e_dashde
