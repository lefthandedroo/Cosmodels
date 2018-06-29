7#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:38:48 2018

@author: BallBlueMeercat
"""
from scipy.integrate import odeint
from firstderivs import edecay, Hdecay, rdecay_m, rdecay_de, interacting, LCDM
import numpy as np

firstderivs_functions = {'edecay':edecay,
                         'Hdecay':Hdecay,
                         'rdecay_m':rdecay_m,
                         'rdecay_de':rdecay_de,
                         'interacting':interacting,
                         'LCDM':LCDM}

def zodesolve(params, zpicks, firstderivs_key):
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
    
    # Standard cosmological parameters.
    H0 = 1
    c_over_H0 = 4167 * 10**6    # c/H0 in parsecs
    
    # Initial conditions at z = 0 (now).
    t0 = 0              # time
    a0 = 1.0            # scale factor
    z0 = 0              # redshift
    dl0 = 0             # luminosity distance
    rho_c0 = H0**2      # critical density
    ombar_m0 = params.get('m', 0)                        # e_m(z)/ec(z=0)
    gamma = params.get('gamma',0)
    ombar_de0 = params.get('de', rho_c0/rho_c0 -ombar_m0)# e_de(z)/ec(z=0)
    
    # ODE solver parameters:
    abserr = 1.0e-8
    relerr = 1.0e-6
    
    # Pack up the initial conditions and eq of state parameters.
    v0 = [t0, a0, ombar_m0, ombar_de0, z0, dl0]
    
    # Extracting the parsed mode of interaction.
    firstderivs_function = firstderivs_functions.get(firstderivs_key,0)
    if firstderivs_function == 0:
        print("firstderivs_functions dict didn't have the key zodeosolve asked for")
    
    # Call the ODE solver. maxstep=5000000 added later to try and avoid
    if firstderivs_key == 'LCDM':
        vsol = odeint(firstderivs_function, v0, zpicks, args=(H0,), 
                  atol=abserr, rtol=relerr, mxstep=5000000)  
    else:
        vsol = odeint(firstderivs_function, v0, zpicks, args=(gamma,H0), 
                  atol=abserr, rtol=relerr, mxstep=5000000)
            
    # Separate results into their own arrays:
    t = vsol[1:,0]
    a = vsol[1:,1]
    ombar_m = vsol[1:,2]
    ombar_de = vsol[1:,3]
    z = vsol[1:,4]    
    dl = vsol[1:,5] * (1+z)  # in units of dl*(H0/c)
    dlpc = dl * c_over_H0    # dl in parsecs (= vsol[dl] * c/H0)
    
#    print('ombar_m max = ',max(ombar_m))
#    print('ombar_m min = ', min(ombar_m))
#    print('ombar_de max = ',max(ombar_de))
#    print('ombar_de min = ', min(ombar_de))
    
    plot_var = t, dlpc, dl, a, ombar_m, gamma, ombar_de, ombar_m0, ombar_de0
    
    return dlpc, plot_var

def odesolve(params, zpicks, firstderivs_key):
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
    
    # Standard cosmological parameters.
    H0 = 1
    c_over_H0 = 4167 * 10**6    # c/H0 in parsecs
    
    # Initial conditions at z = 0 (now).
    dl0 = 0             # luminosity distance
    rho_c0 = H0**2      # critical density
    ombar_m0 = params.get('m', 0)                        # e_m(z)/ec(z=0)
    gamma = params.get('gamma',0)
    ombar_de0 = params.get('de', rho_c0/rho_c0 -ombar_m0)# e_de(z)/ec(z=0)
    
    # ODE solver parameters:
    abserr = 1.0e-8
    relerr = 1.0e-6
    
    # Pack up the initial conditions and eq of state parameters.
    v0 = [ombar_m0, ombar_de0, dl0]
    
    # Extracting the parsed mode of interaction.
    firstderivs_function = firstderivs_functions.get(firstderivs_key,0)
    if firstderivs_function == 0:
        print("firstderivs_functions dict didn't have the key zodeosolve asked for")
    
    # Call the ODE solver. maxstep=5000000 added later to try and avoid
    if firstderivs_key == 'LCDM':
        vsol = odeint(firstderivs_function, v0, zpicks, args=(H0,), 
                  atol=abserr, rtol=relerr, mxstep=5000000)  
    else:
        vsol = odeint(firstderivs_function, v0, zpicks, args=(gamma,H0), 
                  atol=abserr, rtol=relerr, mxstep=5000000)
            
    # Separate results into their own arrays:
    z = np.asarray(zpicks)
    z = z[1:]
    ombar_m = vsol[1:,0]
    ombar_de = vsol[1:,1]
    dl = vsol[1:,2] * (1+z)  # in units of dl*(H0/c)
    dlpc = dl * c_over_H0    # dl in parsecs (= vsol[dl] * c/H0)
    
    
    plot_var = dlpc, dl, ombar_m, gamma, ombar_de, ombar_m0, ombar_de0
    
    return dlpc, plot_var