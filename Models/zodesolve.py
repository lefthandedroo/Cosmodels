7#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:38:48 2018

@author: BallBlueMeercat
"""
from scipy.integrate import odeint
import firstderivs_cython as f
#import firstderivs as f


firstderivs_functions = {
        'exotic':f.exotic,
        'late_intxde':f.late_intxde,
        'heaviside_late_int':f.heaviside_late_int,
        'late_int':f.late_int,
        'expgamma':f.expgamma,
        'txgamma':f.txgamma,
        'zxgamma':f.zxgamma,
        'gamma_over_z':f.gamma_over_z,
        'zxxgamma':f.zxxgamma,
        'gammaxxz':f.gammaxxz,
        'rdecay_m':f.rdecay_m,
        'rdecay_de':f.rdecay_de,
        'rdecay_mxde':f.rdecay_mxde,
        'rdecay':f.rdecay,                         
        'interacting':f.interacting,
        'LCDM':f.LCDM
        }

def zodesolve(params, zpicks, firstderivs_key):
    """
    Takes in:
        params = dictionary w/ key:value
            'gamma':int/float = interaction constant;
            'm':int/float = e_m(t)/ec(t0) at t=t0;
        zpicks = list of redshifts ;
        firstderivs_key = string, name of IVCDM model to use for model mag.
    
    """
    plot_var = {}
    
    # Inserting 0 at the front of redshifts to use initial conditions.
    zpicks = [0.0] + zpicks
    
    # Standard cosmological parameters.
    H0 = 1.0
    c_over_H0 = 4167 * 10**6    # c/H0 in parsecs
    
    # Initial conditions at z = 0 (now).
    t0 = 0.0              # time
    a0 = 1.0              # scale factor
    z0 = 0.0              # redshift
    dl0 = 0.0             # luminosity distance
    rho_c0 = H0**2        # critical density
    ombar_r0 = 0.0        # e_r(z)/ec(z=0)
    for key in params[0]:
        ombar_m0 = params[0].get(key, 0)    # e_m(z)/ec(z=0)
        plot_var['ombar_m0'] = ombar_m0
    ombar_de0 = rho_c0/rho_c0 - ombar_m0    # e_de(z)/ec(z=0)
  
    # Packing up interaction terms:
    int_terms = []
    for i in range(4,len(params)):
        for key in params[i]:
            int_terms.append(params[i][key])
            plot_var[key] = params[i][key]
        
    # Extracting the parsed mode of interaction.
    firstderivs_function = firstderivs_functions.get(firstderivs_key,0)
    if firstderivs_function == 0:
        print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
        print("firstderivs_functions didn't have the key zodeosolve parsed")
        print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
    
    # Pack up the initial conditions and eq of state parameters.
    if firstderivs_key == 'exotic':
        v0 = [t0, a0, ombar_m0, ombar_r0, ombar_de0, z0, dl0]
    else:
        v0 = [t0, a0, ombar_m0, ombar_de0, z0, dl0]
    
    # Call the ODE solver. 
    vsol = odeint(firstderivs_function, v0, zpicks, args=(int_terms,H0), 
                  atol=1.0e-8, rtol=1.0e-6)        

    # Separate results into their own arrays:
    plot_var['t'] = vsol[1:,0]
    plot_var['a'] = vsol[1:,1]
    plot_var['ombar_m'] = vsol[1:,2]
    if firstderivs_key == 'exotic':  
        plot_var['ombar_r'] = vsol[1:,3]
        plot_var['ombar_de'] = vsol[1:,4]
        plot_var['z'] = vsol[1:,5]    
        plot_var['dl'] = vsol[1:,6] * (1+plot_var['z'])  # in units of dl*(H0/c)
    else:
        plot_var['ombar_de'] = vsol[1:,3]
        plot_var['z'] = vsol[1:,4]    
        plot_var['dl'] = vsol[1:,5] * (1+plot_var['z'])  # in units of dl*(H0/c)
    dlpc = plot_var['dl'] * c_over_H0    # dl in parsecs (= vsol[dl] * c/H0)
    plot_var['ombar_r0'] = ombar_r0
    plot_var['ombar_de0'] = ombar_de0     
    
    return dlpc, plot_var