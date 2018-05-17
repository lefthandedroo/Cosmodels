7#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:38:48 2018

@author: BallBlueMeercat
"""
from scipy.integrate import odeint

import zfirstderivs


# Standard cosmological parameters.
H0 = 1
c_over_H0 = 4167 * 10**6    # c/H0 in parsecs

def zodesolve(params, zpicks):
    """
    Takes in:
        gamma = interaction constant;
        m = e_m(t)/ec(t0) at t=t0;
#        de = e_de(t)/ec(t0) at t=t0.
    Returns: 
        z = numpoints number of redshifts zmin<z<zmax;
        dlpc = luminosity distance in pc.
    
    """
#    print('@@ zodesolve has been called')

    # Inserting 0 at the front of redshifts to allow initial conditions.
    zpicks = [0.0] + zpicks
    
    print('zodesolve, type of params',type(params))
    
    # Initial conditions at z = 0 (now).
    t0 = 0              # time
    a0 = 1.0            # scale factor
    z0 = 0              # redshift
    dl0 = 0             # luminosity distance
    rho_c0 = H0**2      # critical density
    ombar_m0 = params.get('m_true', 0)                        # e_m(z)/ec(z=0)
    ombar_de0 = params.get('de_true', rho_c0/rho_c0 -ombar_m0)# e_de(z)/ec(z=0)
    gamma = params.get('gamma_true',0)
    
    # ODE solver parameters:
    abserr = 1.0e-8
    relerr = 1.0e-6
    
    # Pack up the initial conditions and eq of state parameters.
    v0 = [t0, a0, ombar_m0, ombar_de0, z0, dl0]
    
    # Call the ODE solver. maxstep=5000000 added later to try and avoid 
    vsol = odeint(zfirstderivs.zfirstderivs, v0, zpicks, args=(gamma,H0), 
                  atol=abserr, rtol=relerr, mxstep=5000000)
            
    # Separate results into their own arrays:
    t = vsol[1:,0]
    a = vsol[1:,1]
    ombar_m = vsol[1:,2]
    ombar_de = vsol[1:,3]
    z = vsol[1:,4]    
    dl = vsol[1:,5] * (1+z)   # in units of dl*(H0/c)
    dlpc = dl * c_over_H0    # dl in parsecs (= vsol[dl] * c/H0)
    
    return t, dlpc, dl, a, ombar_m, ombar_de, ombar_de0
