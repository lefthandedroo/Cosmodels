7#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:38:48 2018

@author: BallBlueMeercat
"""
import numpy as np
from scipy.integrate import odeint

import zfirstderivs

#from matplotlib.font_manager import FontProperties
#from pylab import figure, plot, xlabel, ylabel, grid, legend, title, show, axis


# Standard cosmological parameters.
H0 = 1       # Hubble parameter at t=now
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
#    print('@@ odesolve has been called')
 
    # Initial conditions at z = 0.
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
    
    # Create z samples for the ODE solver.
    
#    import random
#    zpicks = random.sample(range(0, numpoints*8), numpoints-1)
#    zpicks.sort()
#    denom = numpoints*10
#    zpicks = [time / denom for time in zpicks]
#    zpicks = [0] + zpicks
#    zpicks = [ -x for x in zpicks]
#    print('len zpicks after creation', len(zpicks))
#    print('zpicks:')
#    print(zpicks)
    
    # Pack up the initial conditions and eq of state parameters.
    v0 = [a0, a_dot0, e_dash0m, e_dash0de, z0, dl0]
    
    # Call the ODE solver. maxstep=5000000 added later to try and avoid 
    vsol = odeint(zfirstderivs.zfirstderivs, v0, zpicks, args=(gamma,), 
                  atol=abserr, rtol=relerr, mxstep=5000000)
            
    # Remove unwanted results which are too close to big bang from the plot.
    # Separate results into their own arrays:
    a = vsol[:,0]
    a_dot = vsol[:,1]
    t = zpicks
    t_cut = t
    a_cut = a
    a_dotcut = a_dot
    e_dashm = vsol[:,2]
    e_dashde = vsol[:,3]
    z = vsol[:,4]
    dl = vsol[:,5] * (1+z)   # in units of dl*(H0/c)
    dlpc = dl * c_over_H0    # dl in parsecs (= vsol[dl] * c/H0)

        
    # Take out values at the first redshift z=0:
    blowup = 1
    t_cut = np.asarray(t)
    
    a = a[blowup:]
    a_dot = a_dot[blowup:]
    t = t[blowup:]
    t_cut = t_cut[blowup:]
    a_cut = a_cut[blowup:]
    a_dotcut = a_dotcut[blowup:]
    e_dashm = e_dashm[blowup:]
    e_dashde = e_dashde[blowup:]
    z = z[blowup:]
    dl = dl[blowup:]
    dlpc = dlpc[blowup:]
    zpicks = zpicks[blowup:]

#    # Age of the universe.
#    age = t_cut[np.argmin(t_cut)]
#    age = -round(age, 2)
    
    
    return z, dlpc, dl, gamma, e_dash0m, e_dash0de, t, a, a_dot, t_cut, a_cut, a_dotcut, e_dashm, e_dashde
